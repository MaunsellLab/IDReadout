  function [allProbeDirs, staleProbeDirs] = makeKernels(replace)
% makeKernels  Build all ordinary per-probe derived analysis products.
%
% makeKernels is the authoritative producer of analysis headers:
%
%   sessionHeader
%       Parent-session metadata shared by all probe splits.
%       Includes timing, preferred direction/noise amplitude, and the
%       session-level probe manifest:
%           probeDirectionsDeg
%           probeTags
%
%   sessionProbeHeader
%       Probe-specific metadata for one derived probe session.
%       Includes probeDirDeg, probeTag, trialIdx, probeCohNoisePC, and
%       probe-session trial counts.
%
% Derived per-probe files written by makeKernels contain both headers:
%   Data/probeXX/NoiseMatrices
%   Data/probeXX/Kernels
%   Data/probeXX/KernelSummaries
%
% Aggregate files, such as average kernels and across-offset summaries, do
% not contain a single sessionHeader/sessionProbeHeader because they span
% multiple sessions.

if nargin < 1 || isempty(replace)
  replace = false;
end
path = folderPath();
staleProbeDirs = [];
allProbeDirs = [];

% ---- Define directories and create as needed ----
[dataFolder, existed] = validFolder(fullfile(path, 'Data', 'Converted'));
if ~existed
  error('makeKernels: MissingConvertedFolder -- Converted data folder not found: %s', dataFolder);
end
[plotRoot] = validFolder(fullfile(path, 'Plots', 'Kernels'));

% ---- Find all relevant .mat data files ----
allMatFiles = dir(fullfile(dataFolder, '*.mat'));
if isempty(allMatFiles)
  fprintf('No .mat files found in %s\n', dataFolder);
  return;
end
names = {allMatFiles.name};
isFileInfo = endsWith(names, '_fileInfo.mat');
dataFiles = allMatFiles(~isFileInfo);
if isempty(dataFiles)
  fprintf('No data .mat files (excluding *_fileInfo.mat) found in %s\n', dataFolder);
  return;
end
[~, sideTypeNames] = sideTypeIndex();

% ---- Process each data file ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  dataFilePath = fullfile(dataFolder, dataFileName);
  [~, baseName] = fileparts(dataFileName);
  clear trials;

  % Load lightweight parent metadata first.  sessionHeader is the compact
  % analysis-facing manifest used for stale-output checks.
  load(dataFilePath, 'header', 'sessionHeader');

  needSessionHeaderRefresh = replace || ...
    ~exist('sessionHeader', 'var') || isempty(sessionHeader) || ...
    ~isfield(sessionHeader, 'probeDirectionsDeg') || isempty(sessionHeader.probeDirectionsDeg) || ...
    ~isfield(sessionHeader, 'probeTags') || isempty(sessionHeader.probeTags);

  if needSessionHeaderRefresh
    load(dataFilePath, 'trials');
    trialMeta = trialMetaFromTrials(header, trials);
    sessionHeader = makeSessionHeader(header, trialMeta);
    save(dataFilePath, 'sessionHeader', '-append');
  end

  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;

  % ---- Tally probe directions and check whether any outputs are missing ----
  needsKernels = replace;

  for p = 1:numel(probeDirectionsDeg)
    probeDirDeg = probeDirectionsDeg(p);
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>

    probeTag = char(probeTags{p});
    probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
    probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
    kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
    matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));
    analysisBaseName = sprintf('%s_%s', baseName, probeTag);
    kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
    matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
    plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

    if ~replace
      outputsStale = ~isfile(plotFilePath) || ~isfile(kernelFilePath) || ~isfile(matrixFilePath);
      if outputsStale
        staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
        needsKernels = true;
      end
    end
  end
  if ~needsKernels
    continue
  end
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbeDirection(header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeDirDeg = probeSessions(p).probeDirDeg;
    probeTag = probeSessions(p).probeTag;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionHeader = probeSessions(p).sessionHeader;
    probeTrials = probeSessions(p).trials;
    analysisBaseName = sprintf('%s_%s', baseName, probeTag);

    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
      extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, probeTrials, [1 2]);
    lr = sessionLRMap(probeTrials);

    sessionData = struct;
    sessionData.sessionHeader = sessionHeader;
    sessionData.sessionProbeHeader = sessionProbeHeader;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);

    summary = compileKernelSummary(sessionHeader, sessionProbeHeader, kernels, probeTag);
    
    save(matrixFilePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr','prefNoiseByPatch', 'probeNoiseByPatch', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    save(kernelFilePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    % summaryFolder = validFolder(fullfile(probeDataFolder, 'KernelSummaries'));
    % 
    % compileKernelSummary(kernelFilePath, ...
    %   'noiseFile', matrixFilePath, 'summaryDir', summaryFolder, 'replace', true);

    titleStr = sprintf('\\bf%d° Probe Kernels %s', probeDirDeg, analysisBaseName);
    plotKernels(1, titleStr, sessionHeader, kernels(1:5,:,:,:), kVars(1:5,:,:), compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end

function trialMeta = trialMetaFromTrials(header, trials)

nTrials = numel(trials);
trialProbeDirs = nan(1, nTrials);
trialHasNoise = false(1, nTrials);

for t = 1:nTrials
  if isfield(trials{t}, 'trial') && isfield(trials{t}.trial, 'data')
    D = trials{t}.trial.data;

    if isfield(D, 'cohNoise')
      trialHasNoise(t) = logical(D.cohNoise);
    else
      trialHasNoise(t) = true;
    end

    if isfield(D, 'probeDirDeg')
      trialProbeDirs(t) = double(D.probeDirDeg);
    end
  end
end

% Old single-probe compatibility at the conversion boundary only.
if all(isnan(trialProbeDirs(trialHasNoise)))
  if isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
    trialProbeDirs(:) = double(header.probeDirDeg.data);
  else
    error('trialMetaFromTrials:MissingProbeDir', ...
      'No per-trial probeDirDeg and no header.probeDirDeg.data.');
  end
end

if any(isnan(trialProbeDirs(trialHasNoise)))
  error('trialMetaFromTrials:IncompleteProbeDir', ...
    'Missing probeDirDeg for at least one noise trial.');
end

probeDirectionsDeg = unique(trialProbeDirs(trialHasNoise));
probeDirectionsDeg(probeDirectionsDeg == -1) = [];

trialMeta = struct();
trialMeta.probeDirectionsDeg = probeDirectionsDeg(:)';
trialMeta.nProbeDirections = numel(probeDirectionsDeg);
trialMeta.probeTags = arrayfun(@(d) sprintf('probe%d', round(d)), ...
  trialMeta.probeDirectionsDeg, 'UniformOutput', false);
trialMeta.nNoiseTrials = sum(trialHasNoise & trialProbeDirs ~= -1);
end

%% compileKernelSummary
function summary = compileKernelSummary(sessionHeader, sessionProbeHeader, kernels, probeTag)
% compileKernelSessionSummary
% Build a deterministic lightweight summary for one kernel file.

% P = inputParser;
% addParameter(P, 'noiseFile', '', @ischar);
% addParameter(P, 'summaryDir', '', @ischar);
% addParameter(P, 'replace', false, @islogical);
% addParameter(P, 'trackSideType', 2, @(x) isnumeric(x) && isscalar(x));
% addParameter(P, 'trackStepType', 2, @(x) isnumeric(x) && isscalar(x));
% addParameter(P, 'minTrials', 200, @(x) isnumeric(x) && isscalar(x));
% addParameter(P, 'minPrefEnergy', 1e-6, @(x) isnumeric(x) && isscalar(x));
% parse(P, varargin{:});
% R = P.Results;

% [kDir, baseName, ~] = fileparts(kernelFile);
% 
% if isempty(R.summaryDir)
%   dataDir = fileparts(kDir);
%   summaryDir = fullfile(dataDir, 'KernelSummaries');
% else
%   summaryDir = R.summaryDir;
% end
% if ~exist(summaryDir, 'dir')
%   mkdir(summaryDir);
% end
% 
% summaryFile = fullfile(summaryDir, [baseName '_kernelSummary.mat']);
% 
% if exist(summaryFile, 'file') && ~R.replace
%   tmp = load(summaryFile, 'summary');
%   summary = tmp.summary;
%   return
% end

% K = load(kernelFile);
% sessionHeader = K.sessionHeader;
% sessionProbeHeader = K.sessionProbeHeader;

frameRateHz = localDataValue(sessionHeader.frameRateHz);
msPerVFrame = 1000 / frameRateHz;

[preStepMS, intStartMS, intDurMS] = integralWindowMS();
iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
  round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
iIndices = iIndices(:)';

nTime = size(kernels, 4);
iIndices = iIndices(iIndices >= 1 & iIndices <= nTime);

% if isempty(R.noiseFile)
%   noiseDir = fullfile(fileparts(kDir), 'NoiseMatrices');
%   noiseFile = fullfile(noiseDir, [baseName '.mat']);
% else
%   noiseFile = R.noiseFile;
% end

summary = struct;
summary.version = 3;
[~, baseName] = fileparts(sessionHeader.fileName);
summary.sessionID = baseName;

probeDataFolder = fullfile(folderPath, 'Data', probeTag);
kernelFolder = fullfile(probeDataFolder, 'Kernels');
noiseFolder = fullfile(probeDataFolder, 'NoiseMatrices');

% analysisBaseName = sprintf('%s_%s', sessionHeader.fileName, probeTag);
summary.kernelFile = fullfile(kernelFolder, [baseName, '.mat']);
summary.noiseFile = fullfile(noiseFolder, [baseName, '.mat']);
summary.probeOffsetDeg = sessionProbeHeader.probeDirDeg;

% summary.track = struct;
% summary.track.sideType = R.trackSideType;
% summary.track.stepType = R.trackStepType;
% summary.track.sideLabel = localSideTypeLabel(R.trackSideType);
% summary.track.stepLabel = localStepTypeLabel(R.trackStepType);

% summary.metrics = struct;
% summary.metrics.nTrialsByStep = hitStats.nTrials;
% summary.metrics.nHitsByStep   = hitStats.nHits;
% summary.metrics.nTrials       = sum(hitStats.nTrials);
% summary.metrics.nHits         = sum(hitStats.nHits);
% summary.metrics.pCorrect      = 100 * summary.metrics.nHits / max(summary.metrics.nTrials, 1);
% summary.metrics.nRFTrials     = hitStats.nRFTrials;
% summary.metrics.nRFHits       = hitStats.nRFHits;

% summary.scale = struct;
% % Raw scale: direct probe/pref fit using the measured kernel amplitudes.
% if isfield(compStats, 'rawScale')
%   summary.scale.rawEstimate = compStats.rawScale(R.trackSideType, R.trackStepType);
%   summary.scale.rawSEM      = compStats.rawScaleSEM(R.trackSideType, R.trackStepType);
% % else
% %   % Backward compatibility for older kernel files.
% %   summary.scale.rawEstimate = compStats.scale(R.trackSideType, R.trackStepType);
% %   summary.scale.rawSEM      = compStats.scaleSEM(R.trackSideType, R.trackStepType);
% end

% Normalized scale: probe kernel rescaled to the pref-noise amplitude
% convention before fitting. This is the primary tracked value.
% if isfield(compStats, 'normScale')
%   summary.scale.normEstimate = compStats.normScale(R.trackSideType, R.trackStepType);
%   summary.scale.normSEM      = compStats.normScaleSEM(R.trackSideType, R.trackStepType);
% % else
%   % Older kernel files did not store normalized values. Fall back loudly but
%   % retain compatibility so historical files can still be inspected.
%   warning('compileKernelSessionSummary:MissingNormScale', ...
%     'Kernel file %s lacks compStats.normScale; using raw scale as normalized fallback.', baseName);
%   summary.scale.normEstimate = summary.scale.rawEstimate;
%   summary.scale.normSEM      = summary.scale.rawSEM;
% end

% summary.integral = struct;
% summary.ratio = struct;
% 
% if isfield(compStats, 'rawIntegrals')
%   summary.integral.raw = squeeze(compStats.rawIntegrals(R.trackSideType, R.trackStepType, :))';
%   summary.ratio.raw    = compStats.rawR(R.trackSideType, R.trackStepType);
% else
%   summary.integral.raw = squeeze(compStats.kIntegrals(R.trackSideType, R.trackStepType, :))';
%   summary.ratio.raw    = compStats.R(R.trackSideType, R.trackStepType);
% end
% 
% if isfield(compStats, 'normIntegrals')
%   summary.integral.norm = squeeze(compStats.normIntegrals(R.trackSideType, R.trackStepType, :))';
%   summary.ratio.norm    = compStats.normR(R.trackSideType, R.trackStepType);
% else
%   summary.integral.norm = summary.integral.raw;
%   summary.ratio.norm    = summary.ratio.raw;
% end
% 
% if isfield(compStats, 'normInfo')
%   summary.normInfo = compStats.normInfo;
% else
%   summary.normInfo = struct();
% end

% Primary aliases used by downstream code.
% summary.scale.estimate = summary.scale.normEstimate;
% summary.scale.sem      = summary.scale.normSEM;
% 
% summary.scale.valid      = false;
% kPref = squeeze(kernels(R.trackSideType, R.trackStepType, 1, iIndices));
% summary.pref = struct;
% summary.pref.energy   = sum(kPref(:).^2);
% summary.pref.integral = sum(kPref(:)) * msPerVFrame;
% summary.pref.nBins    = numel(iIndices);
% % 
% summary.bootstrap = struct;
% summary.bootstrap.done    = false;
% summary.bootstrap.nReps   = 0;
% summary.bootstrap.method  = 'trial';
% summary.bootstrap.rngSeed = NaN;

summary.flags = struct;
% summary.flags.lowTrialCount    = sum(hitStats.nTrials) < R.minTrials;
% summary.flags.lowPrefEnergy    = summary.pref.energy < R.minPrefEnergy;
% summary.flags.unstableScale    = false;
% summary.flags.missingNoiseFile = ~exist(noiseFile, 'file');
% summary.flags.needsRefresh     = false;

% summary.scale.valid = ...
%   isfinite(summary.scale.estimate) && ...
%   ~summary.flags.lowTrialCount && ...
%   ~summary.flags.lowPrefEnergy && ...
%   ~summary.flags.unstableScale;

end

% function dt = localDateFromHeaderOrFile(header, baseName)
% dt = NaT;
% 
% try
%   if isfield(header, 'date')
%     raw = header.date.data;
%     if isnumeric(raw)
%       dt = datetime(raw, 'ConvertFrom', 'datenum');
%       return
%     elseif ischar(raw) || isstring(raw)
%       dt = datetime(raw);
%       return
%     end
%   end
% catch
% end
% 
% tok = regexp(baseName, '\d{8}', 'match', 'once');
% if ~isempty(tok)
%   dt = datetime(tok, 'InputFormat', 'yyyyMMdd');
% end
% end

% function txt = localStepTypeLabel(stepType)
% switch stepType
%   case 1
%     txt = 'DEC';
%   case 2
%     txt = 'INC';
%   otherwise
%     txt = sprintf('step %d', stepType);
% end
% end

% function txt = localSideTypeLabel(sideType)
% switch sideType
%   case 1
%     txt = 'change-noChange';
%   case 2
%     txt = 'changeSide';
%   case 3
%     txt = 'noChangeSide';
%   case 4
%     txt = 'RF';
%   case 5
%     txt = 'Opp';
%   otherwise
%     txt = sprintf('side %d', sideType);
% end
% end

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