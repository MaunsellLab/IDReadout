function summary = compileKernelSessionSummary(kernelFile, varargin)
% compileKernelSessionSummary
%
% Build and optionally save a lightweight per-session summary for tracking
% one selected kernel scale across sessions.
%
% Uses:
%   - Kernels file for observed metrics
%   - NoiseMatrices file for within-session trial bootstrap
%
% Default tracked entry:
%   sideType = 2   ('changeSide')
%   stepType = 2   ('INC')
%
% Example:
%   compileKernelSessionSummary( ...
%     '/Users/Shared/Data/IDReadout/Data/Kernels/IDReadout_Meetz_20260402.mat', ...
%     'doBootstrap', true, 'nBoot', 500);

P = inputParser;
addParameter(P, 'noiseFile', '', @ischar);
addParameter(P, 'summaryDir', '', @ischar);
addParameter(P, 'replace', false, @islogical);
addParameter(P, 'doBootstrap', true, @islogical);
addParameter(P, 'nBoot', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(P, 'trackSideType', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'trackStepType', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'minTrials', 200, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'minPrefEnergy', 1e-6, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'rngSeed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(P, varargin{:});
R = P.Results;

[kDir, baseName, ~] = fileparts(kernelFile);

if isempty(R.summaryDir)
  dataDir = fileparts(kDir);
  summaryDir = fullfile(dataDir, 'KernelSummaries');
else
  summaryDir = R.summaryDir;
end
if ~exist(summaryDir, 'dir')
  mkdir(summaryDir);
end

summaryFile = fullfile(summaryDir, [baseName '_kernelSummary.mat']);

if exist(summaryFile, 'file') && ~R.replace
  tmp = load(summaryFile, 'summary');
  summary = tmp.summary;
  return
end

K = load(kernelFile);

frameRateHz = K.header.frameRateHz.data;
msPerVFrame = 1000 / frameRateHz;

[preStepMS, intStartMS, intDurMS] = integralWindowMS();
iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
  round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
iIndices = iIndices(:)';

nTime = size(K.kernels, 4);
iIndices = iIndices(iIndices >= 1 & iIndices <= nTime);

if isempty(R.noiseFile)
  noiseDir = fullfile(fileparts(kDir), 'NoiseMatrices');
  noiseFile = fullfile(noiseDir, [baseName '.mat']);
else
  noiseFile = R.noiseFile;
end

summary = struct;
summary.version = 1;
summary.sessionID = baseName;
summary.date = localDateFromHeaderOrFile(K.header, baseName);
summary.kernelFile = kernelFile;
summary.noiseFile = noiseFile;

summary.track = struct;
summary.track.sideType = R.trackSideType;
summary.track.stepType = R.trackStepType;
summary.track.sideLabel = localSideTypeLabel(R.trackSideType);
summary.track.stepLabel = localStepTypeLabel(R.trackStepType);

summary.metrics = struct;
summary.metrics.nTrialsByStep = K.hitStats.nTrials;
summary.metrics.nHitsByStep   = K.hitStats.nHits;
summary.metrics.nTrials       = sum(K.hitStats.nTrials);
summary.metrics.nHits         = sum(K.hitStats.nHits);
summary.metrics.pCorrect      = 100 * summary.metrics.nHits / max(summary.metrics.nTrials, 1);
summary.metrics.nRFTrials     = K.hitStats.nRFTrials;
summary.metrics.nRFHits       = K.hitStats.nRFHits;

summary.scale = struct;
summary.scale.estimate   = K.compStats.scale(R.trackSideType, R.trackStepType);
summary.scale.sem        = K.compStats.scaleSEM(R.trackSideType, R.trackStepType);
summary.scale.ci68       = [NaN NaN];
summary.scale.ci95       = [NaN NaN];
summary.scale.ci68Width  = NaN;
summary.scale.ci95Width  = NaN;
summary.scale.bootMedian = NaN;
summary.scale.bootSD     = NaN;
summary.scale.valid      = false;

kPref = squeeze(K.kernels(R.trackSideType, R.trackStepType, 1, iIndices));
summary.pref = struct;
summary.pref.energy   = sum(kPref(:).^2);
summary.pref.integral = sum(kPref(:)) * msPerVFrame;
summary.pref.nBins    = numel(iIndices);

summary.bootstrap = struct;
summary.bootstrap.done    = false;
summary.bootstrap.nReps   = 0;
summary.bootstrap.method  = 'trial';
summary.bootstrap.rngSeed = NaN;

summary.flags = struct;
summary.flags.excluded         = excludeFile(K.header);
summary.flags.lowTrialCount    = summary.metrics.nTrials < R.minTrials;
summary.flags.lowPrefEnergy    = summary.pref.energy < R.minPrefEnergy;
summary.flags.unstableScale    = false;
summary.flags.missingNoiseFile = ~exist(noiseFile, 'file');
summary.flags.needsRefresh     = false;

if R.doBootstrap && exist(noiseFile, 'file')

  if ~isempty(R.rngSeed)
    rng(R.rngSeed);
    summary.bootstrap.rngSeed = R.rngSeed;
  end

  N = load(noiseFile);

  try
    prefCohNoisePC  = K.header.prefNoiseCohPC.data;
    probeCohNoisePC = K.header.probeNoiseCohPC.data;
    % [prefCohNoisePC, probeCohNoisePC] = localGetNoiseAmplitudes(K.header);
  catch ME
    warning('compileKernelSessionSummary:noiseAmplitudes', ...
      'Could not extract noise amplitudes from header for %s: %s', ...
      baseName, ME.message);
    prefCohNoisePC = NaN;
    probeCohNoisePC = NaN;
  end

  nTrials = size(N.prefNoiseByPatch, 3);
  bootScale = nan(R.nBoot, 1);

  for b = 1:R.nBoot
    idx = randi(nTrials, nTrials, 1);

    try
      [kernelsBoot, ~, ~] = computeSessionKernelsFromNoiseMatrices( ...
        N.prefNoiseByPatch(:,:,idx), ...
        N.probeNoiseByPatch(:,:,idx), ...
        N.trialOutcomesAll(idx), ...
        N.changeSidesAll(idx), ...
        N.changeIndicesAll(idx), ...
        prefCohNoisePC, probeCohNoisePC);

      [scaleBoot, ~, ~, ~] = kernelScaleFit(kernelsBoot, msPerVFrame);
      bootScale(b) = scaleBoot(R.trackSideType, R.trackStepType);

    catch
      bootScale(b) = NaN;
    end
  end

  bootScale = bootScale(isfinite(bootScale));

  if ~isempty(bootScale)
    summary.scale.ci68       = prctile(bootScale, [16 84]);
    summary.scale.ci95       = prctile(bootScale, [2.5 97.5]);
    summary.scale.ci68Width  = diff(summary.scale.ci68);
    summary.scale.ci95Width  = diff(summary.scale.ci95);
    summary.scale.bootMedian = median(bootScale);
    summary.scale.bootSD     = std(bootScale, 0);

    summary.bootstrap.done  = true;
    summary.bootstrap.nReps = numel(bootScale);
  else
    summary.flags.unstableScale = true;
  end
end

summary.scale.valid = ...
  isfinite(summary.scale.estimate) && ...
  ~summary.flags.lowTrialCount && ...
  ~summary.flags.lowPrefEnergy && ...
  ~summary.flags.unstableScale;

save(summaryFile, 'summary');

end


function dt = localDateFromHeaderOrFile(header, baseName)
dt = NaT;

try
  if isfield(header, 'date')
    raw = header.date.data;
    if isnumeric(raw)
      dt = datetime(raw, 'ConvertFrom', 'datenum');
      return
    elseif ischar(raw) || isstring(raw)
      dt = datetime(raw);
      return
    end
  end
catch
end

tok = regexp(baseName, '\d{8}', 'match', 'once');
if ~isempty(tok)
  dt = datetime(tok, 'InputFormat', 'yyyyMMdd');
end
end


function txt = localStepTypeLabel(stepType)
switch stepType
  case 1
    txt = 'DEC';
  case 2
    txt = 'INC';
  otherwise
    txt = sprintf('step %d', stepType);
end
end


function txt = localSideTypeLabel(sideType)
switch sideType
  case 1
    txt = 'change-noChange';
  case 2
    txt = 'changeSide';
  case 3
    txt = 'noChangeSide';
  case 4
    txt = 'RF';
  case 5
    txt = 'Opp';
  otherwise
    txt = sprintf('side %d', sideType);
end
end