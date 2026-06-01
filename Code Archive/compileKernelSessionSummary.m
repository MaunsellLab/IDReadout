function summary = compileKernelSessionSummary(kernelFile, varargin)
% compileKernelSessionSummary
%
% Build a deterministic lightweight summary from one kernel file.
% The output file contains:
%   summary
%   sessionHeader
%   sessionProbeHeader
%
% KernelSummary files do not store per-session bootstrap results. Bootstrap
% uncertainty for the readout analysis is computed fresh by
% updateAcrossOffsetSummaries using hierarchical resampling.

P = inputParser;
addParameter(P, 'noiseFile', '', @ischar);
addParameter(P, 'summaryDir', '', @ischar);
addParameter(P, 'replace', false, @islogical);
addParameter(P, 'trackSideType', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'trackStepType', 2, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'minTrials', 200, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'minPrefEnergy', 1e-6, @(x) isnumeric(x) && isscalar(x));
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

sessionHeader = K.sessionHeader;
sessionProbeHeader = K.sessionProbeHeader;

frameRateHz = localDataValue(sessionHeader.frameRateHz);
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
summary.version = 2;
summary.sessionID = baseName;
summary.date = localDateFromHeaderOrFile(sessionHeader, baseName);
summary.kernelFile = kernelFile;
summary.noiseFile = noiseFile;
summary.probeOffsetDeg = sessionProbeHeader.probeDirDeg;

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
% Raw scale: direct probe/pref fit using the measured kernel amplitudes.
if isfield(K.compStats, 'rawScale')
  summary.scale.rawEstimate = K.compStats.rawScale(R.trackSideType, R.trackStepType);
  summary.scale.rawSEM      = K.compStats.rawScaleSEM(R.trackSideType, R.trackStepType);
else
  % Backward compatibility for older kernel files.
  summary.scale.rawEstimate = K.compStats.scale(R.trackSideType, R.trackStepType);
  summary.scale.rawSEM      = K.compStats.scaleSEM(R.trackSideType, R.trackStepType);
end

% Normalized scale: probe kernel rescaled to the pref-noise amplitude
% convention before fitting. This is the primary tracked value.
if isfield(K.compStats, 'normScale')
  summary.scale.normEstimate = K.compStats.normScale(R.trackSideType, R.trackStepType);
  summary.scale.normSEM      = K.compStats.normScaleSEM(R.trackSideType, R.trackStepType);
else
  % Older kernel files did not store normalized values. Fall back loudly but
  % retain compatibility so historical files can still be inspected.
  warning('compileKernelSessionSummary:MissingNormScale', ...
    'Kernel file %s lacks compStats.normScale; using raw scale as normalized fallback.', baseName);
  summary.scale.normEstimate = summary.scale.rawEstimate;
  summary.scale.normSEM      = summary.scale.rawSEM;
end

summary.integral = struct;
summary.ratio = struct;

if isfield(K.compStats, 'rawIntegrals')
  summary.integral.raw = squeeze(K.compStats.rawIntegrals(R.trackSideType, R.trackStepType, :))';
  summary.ratio.raw    = K.compStats.rawR(R.trackSideType, R.trackStepType);
else
  summary.integral.raw = squeeze(K.compStats.kIntegrals(R.trackSideType, R.trackStepType, :))';
  summary.ratio.raw    = K.compStats.R(R.trackSideType, R.trackStepType);
end

if isfield(K.compStats, 'normIntegrals')
  summary.integral.norm = squeeze(K.compStats.normIntegrals(R.trackSideType, R.trackStepType, :))';
  summary.ratio.norm    = K.compStats.normR(R.trackSideType, R.trackStepType);
else
  summary.integral.norm = summary.integral.raw;
  summary.ratio.norm    = summary.ratio.raw;
end

if isfield(K.compStats, 'normInfo')
  summary.normInfo = K.compStats.normInfo;
else
  summary.normInfo = struct();
end

% Primary aliases used by downstream code.
summary.scale.estimate = summary.scale.normEstimate;
summary.scale.sem      = summary.scale.normSEM;


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
summary.flags.lowTrialCount    = summary.metrics.nTrials < R.minTrials;
summary.flags.lowPrefEnergy    = summary.pref.energy < R.minPrefEnergy;
summary.flags.unstableScale    = false;
summary.flags.missingNoiseFile = ~exist(noiseFile, 'file');
summary.flags.needsRefresh     = false;

summary.scale.valid = ...
  isfinite(summary.scale.estimate) && ...
  ~summary.flags.lowTrialCount && ...
  ~summary.flags.lowPrefEnergy && ...
  ~summary.flags.unstableScale;

save(summaryFile, 'summary', 'sessionHeader', 'sessionProbeHeader');
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

function v = localDataValue(x)
if isstruct(x) && isfield(x, 'data')
  v = x.data;
else
  v = x;
end

if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end