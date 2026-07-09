function acrossOffsetSummary = updateAcrossOffsetSummaries(dataDir, varargin)
% This function constructs explicit per-session records from those objects,
% applies readout-fit eligibility rules, performs hierarchical bootstrap
% resampling from the source noise-matrix files, and fits the MT-readout
% model across probe offsets.
%   The required data are in the kernel files, which are separated by probe
%   direction.  
% Kernel file contains:
%   summary
%   sessionHeader
%   sessionProbeHeader% 
% Each kernel file has an associated noise file that provide the
% data needed to do bootstraps in computeSessionKernels:
%       sessionProbeHeader
%       prefNoiseByPatch
%       probeNoiseByPatch
%       trialOutcomesAll
%       changeSidesAll
%       changeIndicesAll
%       lr
%       sideTypeNames
%
% -------------------------------------------------------------------------

if nargin < 1 || isempty(dataDir)
  dataDir = fullfile(domainFolder(mfilename('fullpath')), 'Data');
end
opts = parseInputs(dataDir, varargin{:});

% get info about all kernel files meeting criteria
files = {};
fileInfo = table();
for iDir = 1:numel(opts.kernelsDir)
  thisDir = opts.kernelsDir{iDir};
  [theFiles, theFileInfo] = selectAnalysisFiles(thisDir, 'Animal', opts.Animal, 'FilePattern', opts.FilePattern);
  files = [files; theFiles]; %#ok<AGROW>
  fileInfo = [fileInfo; theFileInfo]; %#ok<AGROW>
end
% Make sure there are no duplicates
[files, ia] = unique(files, 'stable');
fileInfo = fileInfo(ia, :);

% get a list of all session summaries
[sessionList] = loadSessionSummaries(files, fileInfo, opts);
% create an empty acrossOffsetSummary
acrossOffsetSummary = initializeAcrossOffsetSummary(opts, sessionList);
acrossOffsetSummary.fileInfo = fileInfo;
if isempty(sessionList)
    warning('No usable session summaries found in %s.', dataDir);
    saveAcrossOffsetSummary(opts, acrossOffsetSummary);
    return;
end
% extract the data for each probe offset
[offsetData, offsetKeys, usedFiles] = buildOffsetData(sessionList, opts);
if isempty(offsetData) || all([offsetData.nSessions] == 0)
    warning('No included sessions remain after exclusions. Returning summary without bootstrap/model fits.');
    acrossOffsetSummary.empirical = struct([]);
    acrossOffsetSummary.bootstrap = struct();
    acrossOffsetSummary.modelFits = struct();
    acrossOffsetSummary.history = struct([]);
    saveAcrossOffsetSummary(opts, acrossOffsetSummary);
    return;
end
acrossOffsetSummary.type = 'Kernels';
acrossOffsetSummary.offsetData = offsetData;
acrossOffsetSummary.meta.offsetKeysDeg = offsetKeys;
acrossOffsetSummary.meta.summaryFilesUsed = usedFiles;
% load the empirical data from the summaries, including the ci
empirical = computeEmpiricalSummaries(offsetData, opts);
acrossOffsetSummary.empirical = empirical;

bootstrap = runHierarchicalBootstrap(offsetData, opts);
acrossOffsetSummary.bootstrap = bootstrap;

%% ---- Effective MT-to-choice weighting fit ----
offsetsDegAll   = [empirical.probeOffsetDeg];
pooledScaleAll  = [empirical.pooledScale];
bootstrapVarAll = bootstrap.fitBootstrap.offsetFitVar;
nSessionsAll    = [empirical.nSessions];
assert(numel(bootstrapVarAll) == numel(offsetsDegAll), ...
    'bootstrap offset variances do not match empirical offsets');

readoutFitSummary = fitAcrossOffsetReadout(pooledScaleAll, bootstrapVarAll, offsetsDegAll, ...
    'NSessions', nSessionsAll, ...
    'Bounds', opts.Bounds, ...
    'SigmaMTDeg', 37.5, ...
    'PhiDeg', -180:1:179, ...
    'SourceMeasureType', 'kernelScale', ...
    'SourceSideType', opts.ScaleSideType, ...
    'SourceStepType', opts.ScaleStepType, ...
    'SourceMode', opts.ScaleMode);

% Preserve the historical output fields for downstream compatibility.
acrossOffsetSummary.measurements = struct();
acrossOffsetSummary.measurements.offsetsDeg   = offsetsDegAll(:)';
acrossOffsetSummary.measurements.pooledScale  = pooledScaleAll(:)';
acrossOffsetSummary.measurements.bootstrapVar = bootstrapVarAll(:)';
acrossOffsetSummary.readoutModels = readoutFitSummary.readoutModels;
acrossOffsetSummary.readoutModelComparison = readoutFitSummary.readoutModelComparison;
acrossOffsetSummary.readoutModel = readoutFitSummary.readoutModel;
acrossOffsetSummary.history = updateHistory(acrossOffsetSummary, opts);

saveAcrossOffsetSummary(opts, acrossOffsetSummary);
makeAcrossOffsetPlots(acrossOffsetSummary, opts);

if opts.Verbose
  fprintf('updateAcrossOffsetSummaries: done.\n');
end

end

% ========================================================================
function opts = parseInputs(dataDir, varargin)

% locate all the potential summary directories
d = char(cellstr(dataDir));
kernelDirs = {};
probeDirs = dir(fullfile(d, 'probe*'));
probeDirs = probeDirs([probeDirs.isdir]);
for p = 1:numel(probeDirs)
  candidate = fullfile(probeDirs(p).folder, probeDirs(p).name, 'Kernels');
  if exist(candidate, 'dir')
    kernelDirs{end+1} = candidate; %#ok<AGROW>
  end
end
kernelDirs = unique(kernelDirs, 'stable');

p = inputParser;
p.FunctionName = mfilename;

addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addRequired(p, 'dataDir', @(x) ischar(x) || isstring(x) || iscell(x));
addParameter(p, 'PlotDir',  ...
          fullfile(domainFolder(mfilename('fullpath')), 'Plots', 'AcrossProbes', 'ReadoutFits', 'Kernels'), ...
          @(x) ischar(x) || isstring(x));
addParameter(p, 'NBoot', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CILevels', [68 95], @(x) isnumeric(x) && isvector(x) && all(x > 0) && all(x < 100));
addParameter(p, 'Bin180Into179', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'OffsetField', 'probeOffsetDeg', @(x) ischar(x) || isstring(x));
addParameter(p, 'SessionNameField', 'sessionName', @(x) ischar(x) || isstring(x));
addParameter(p, 'ScaleMode', 'scaleFit', @(x) ischar(x) || isstring(x));
addParameter(p, 'ScaleSideType', 'change', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'ScaleStepType', 'inc', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'Bounds', struct(), @(x) isstruct(x));
addParameter(p, 'RandomSeed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));

parse(p, dataDir, varargin{:});
opts = p.Results;

opts.kernelsDir = kernelDirs;
opts.PlotDir  = char(opts.PlotDir);
opts.FilePattern = char(opts.FilePattern);
opts.OffsetField = char(opts.OffsetField);
opts.SessionNameField = char(opts.SessionNameField);
opts.ScaleMode = char(opts.ScaleMode);

saveFile = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', ...
  sprintf('KernelSummary_%s.mat', opts.Animal));
opts.SaveFile = char(saveFile);
if ~exist(fileparts(opts.SaveFile), 'dir')
  mkdir(fileparts(opts.SaveFile));
end
if opts.MakePlots && ~exist(opts.PlotDir, 'dir')
  mkdir(opts.PlotDir);
end

if ~isempty(opts.RandomSeed)
  rng(opts.RandomSeed);
end
end

% ========================================================================
function [sessionList] = loadSessionSummaries(files, fileInfo, opts)
% Load per-session summaries and apply exclusion logic.

if opts.Verbose
    fprintf('  Found %d selected session summary files across %d folders.\n', ...
        numel(files), numel(opts.kernelsDir));
end
sessionList = repmat(emptySessionStruct(), 0, 1);
for iFile = 1:numel(files)
    thisFile = files{iFile};

    load(thisFile, 'sessionHeader', 'sessionProbeHeader', 'compStats', 'hitStats');
    sessionRecord = emptySessionStruct();
    sessionRecord.sessionHeader = sessionHeader;
    sessionRecord.sessionProbeHeader = sessionProbeHeader;
    sessionRecord.sessionName = sessionHeader.fileName;
    sessionRecord.compStats = compStats;
    sessionRecord.hitStats  = hitStats;
    sessionRecord.sourceFile = thisFile;
    sessionRecord.scalePointEstimate = selectCompStatsEntry(compStats.normScale, opts.ScaleSideType, opts.ScaleStepType);
    sessionRecord.nTrials = sum(hitStats.nTrials);
    sessionRecord.nTrialsByStep = hitStats.nTrials;
    sessionRecord.noiseFile  = sessionProbeHeader.probeSessionPath;

    % If we are binning 180 into 179, make that adjustment in the record
    if opts.Bin180Into179 && sessionProbeHeader.probeDirDeg == 180
      sessionRecord.probeOffsetDeg = 179;
    else
      sessionRecord.probeOffsetDeg = sessionProbeHeader.probeDirDeg;
    end
    if ~isempty(fileInfo) && height(fileInfo) >= iFile
        sessionRecord.fileInfo = fileInfo(iFile, :);
    end

    [tfExclude, reasonStr] = excludeFromReadoutFit(sessionRecord.probeOffsetDeg, ...
              sessionRecord.sessionProbeHeader.probeCohNoisePC, opts.Bin180Into179);
    sessionRecord.isExcluded = tfExclude;
    sessionRecord.excludeReason = reasonStr;
    sessionList(end+1, 1) = sessionRecord; %#ok<AGROW>
end
if opts.Verbose
  nExcluded = sum([sessionList.isExcluded]);
  fprintf('  Loaded %d session summaries (%d excluded, %d included).\n', ...
    numel(sessionList), nExcluded, numel(sessionList) - nExcluded);
  if ~isempty(sessionList)
    fprintf('  Session diagnostics:\n');
    for iS = 1:numel(sessionList)
      incTxt = 'In';
      if sessionList(iS).isExcluded
        incTxt = 'Out';
      end
      fprintf('    %-3s  %-36s  offset=%6.1f  nTrials=%5g', ...
        incTxt, sessionList(iS).sessionName, sessionList(iS).probeOffsetDeg, sessionList(iS).nTrials);
      if sessionList(iS).isExcluded && ~isempty(sessionList(iS).excludeReason)
        fprintf('  reason = %s', sessionList(iS).excludeReason);
      end
      fprintf('\n');
    end
  end
end

end

% ========================================================================
function s = emptySessionStruct()

s = struct( ...
    'probeOffsetDeg', NaN, ...
    'sessionName', '', ...
    'sessionDate', NaT, ...
    'summary', struct(), ...
    'sessionHeader', struct(), ...
    'sessionProbeHeader', struct(), ...
    'noiseFile', '', ...
    'compStats', struct(), ...
    'hitStats', struct(), ...
    'rawSummary', struct(), ...
    'scalePointEstimate', NaN, ...
    'nTrials', NaN, ...
    'nTrialsByStep', NaN, ...
    'isExcluded', false, ...
    'excludeReason', '', ...
    'sourceFile', '', ...
    'fileInfo', table() );
end

% ========================================================================
function [exclude, reason] = excludeFromReadoutFit(probeDirDeg, probeCohNoisePC, Bin180Into179)
%excludeFromReadoutFit  Exclude sessions not eligible for paired-probe readout fitting.
%
% These exclusions are specific to the across-offset readout model. They do
% not imply that the session should be excluded from kernel generation,
% kernel averaging, or other preprocessing summaries.
%
% Inputs are explicit to avoid searching through generic session/header
% structs. The caller is responsible for passing canonical values.

exclude = false;
reason = '';

% Paired-probe offsets are strictly between 0 and 180 deg.
% Exact 0 and exact 180 are single-stream probes.
if Bin180Into179
  limitDeg = 180;
else 
  limitDeg = 179;
end
if ~(isfinite(probeDirDeg) && probeDirDeg > 1 && probeDirDeg <= limitDeg)
    exclude = true;
    reason = sprintf('probeDirDeg %.6g is not a paired-probe offset', probeDirDeg);
    return;
end

% Accept known paired-probe amplitudes. This is intentionally a readout-fit
% eligibility rule, not a global file exclusion rule.
targetProbeNoisePCs = [1, 7, 7.07, 10 / sqrt(2), 10];
tol = 2e-3;

if all(abs(probeCohNoisePC - targetProbeNoisePCs) > tol)
    exclude = true;
    reason = sprintf('probeCohNoisePC %.6g is not paired-probe amplitude', probeCohNoisePC);
    return;
end
end

% ========================================================================
function acrossOffsetSummary = initializeAcrossOffsetSummary(opts, sessionList)

acrossOffsetSummary = struct();
acrossOffsetSummary.meta = struct( ...
    'analysisDate', datetime('now'), ...
    'summaryDirs', {opts.kernelsDir}, ...
    'summaryFilesUsed', {{}}, ...
    'nCandidateFiles', numel(sessionList), ...
    'scaleMetric', opts.ScaleMode, ...
    'ciType', 'bootstrap_percentile', ...
    'ciLevels', opts.CILevels(:)', ...
    'nBoot', opts.NBoot, ...
    'bootstrapType', 'hierarchical_session_trial', ...
    'angleUnits', 'deg', ...
    'normalization', 'probe_kernel_amplitude_normalized_to_pref_noise_convention', ...
    'offsetKeysDeg', [], ...
    'readoutBaselineOption', 'inactive_not_identifiable', ...
    'notes', '' );

acrossOffsetSummary.offsetData = struct([]);
acrossOffsetSummary.empirical  = struct([]);
acrossOffsetSummary.bootstrap  = struct();
acrossOffsetSummary.history    = struct([]);
acrossOffsetSummary.modelFits = struct();
acrossOffsetSummary.modelFitsNote = 'Primary interpretation uses acrossOffsetSummary.readoutModel.';
end

% ========================================================================
function [offsetData, offsetKeys, usedFiles] = buildOffsetData(sessionList, opts)
% Group sessions by probe offset and retain both included and excluded items.

allOffsets = [sessionList.probeOffsetDeg];
allOffsets = allOffsets(isfinite(allOffsets));
offsetKeys = unique(allOffsets(:)');
offsetData = repmat(emptyOffsetStruct(), numel(offsetKeys), 1);
for k = 1:numel(offsetKeys)
  thisOffset = offsetKeys(k);
  idx = [sessionList.probeOffsetDeg] == thisOffset;
  these = sessionList(idx);
  includeMask = ~[these.isExcluded];
  included = these(includeMask);

  offsetData(k).probeOffsetDeg = thisOffset;
  offsetData(k).sessionNames = {these.sessionName};
  offsetData(k).sessionDates = [these.sessionDate];
  offsetData(k).nSessions = sum(includeMask);
  offsetData(k).nSessionsTotal = numel(these);
  offsetData(k).nTrialsBySession = [included.nTrials];
  offsetData(k).nTrialsByStep = sum(vertcat(included.nTrialsByStep), 1, 'omitnan');
  offsetData(k).nTrials = sum(offsetData(k).nTrialsBySession, 'omitnan');
  offsetData(k).scaleBySession = [included.scalePointEstimate];
  offsetData(k).scaleSEMBySession = nan(size(offsetData(k).scaleBySession));
  offsetData(k).scaleCILoBySession = nan(size(offsetData(k).scaleBySession));
  offsetData(k).scaleCIHiBySession = nan(size(offsetData(k).scaleBySession));
  offsetData(k).prefEnergyBySession = nan(size(offsetData(k).scaleBySession));
  offsetData(k).fitR2BySession = nan(size(offsetData(k).scaleBySession));
  offsetData(k).includeMask = includeMask;
  offsetData(k).excludeReasons = {these.excludeReason};
  offsetData(k).sourceSummaryFiles = {these.sourceFile};
  offsetData(k).sessionStructs = included;
end

% if we are binning 179° with 180°, combine them into 179°
% if opts.Bin180Into179
%   index179 = find(offsetKeys == 179, 1);
%   index180 = find(offsetKeys == 180, 1);
%   if ~isempty(index179) && ~isempty(index180)
%     offsetData(index179).sessionNames = [offsetData(index179).sessionNames, offsetData(index180).sessionNames];
%     offsetData(index179).sessionDates = [offsetData(index179).sessionDates, offsetData(index180).sessionDates];
%     offsetData(index179).nSessions = offsetData(index179).nSessions + offsetData(index180).nSessions;
%     offsetData(index179).nSessionsTotal = offsetData(index179).nSessionsTotal + offsetData(index180).nSessionsTotal;
%     offsetData(index179).nTrialsBySession = [offsetData(index179).nTrialsBySession, offsetData(index180).nTrialsBySession];
%     offsetData(index179).nTrialsByStep = ...
%           offsetData(index179).nTrialsByStep + offsetData(index180).nTrialsByStep;
%     offsetData(index179).nTrials = offsetData(index179).nSessionsTotal + offsetData(index180).nTrials;
%     offsetData(index179).scaleBySession = [offsetData(index179).scaleBySession, offsetData(index180).scaleBySession];
% 
%     offsetData(index179).scaleSEMBySession = nan(size(offsetData(index179).scaleBySession));
%     offsetData(index179).scaleSEMBySession = nan(size(offsetData(index179).scaleBySession));
%     offsetData(index179).scaleCILoBySession = nan(size(offsetData(index179).scaleBySession));
%     offsetData(index179).scaleCIHiBySession = nan(size(offsetData(index179).scaleBySession));
%     offsetData(index179).prefEnergyBySession = nan(size(offsetData(index179).scaleBySession));
%     offsetData(index179).fitR2BySession = nan(size(offsetData(index179).scaleBySession));
% 
%     offsetData(index179).includeMask = [offsetData(index179).includeMask, offsetData(index180).includeMask];
%     offsetData(index179).excludeReasons = [offsetData(index179).excludeReasons, offsetData(index180).excludeReasons];
%     offsetData(index179).sourceSummaryFiles = [offsetData(index179).sourceSummaryFiles, offsetData(index180).sourceSummaryFiles];
%     offsetData(index179).sessionStructs = [offsetData(index179).sessionStructs; offsetData(index180).sessionStructs];
% 
%     offsetData(index180) = [];
%     offsetKeys(index180) = [];
%   end
% end

% update summary file list with included files only
includedMask = ~[sessionList.isExcluded];
usedFiles = {sessionList(includedMask).sourceFile};
if ~isempty(usedFiles)
  usedFiles = unique(usedFiles, 'stable');
end
end

% ========================================================================
function s = emptyOffsetStruct()

s = struct( ...
    'probeOffsetDeg', NaN, ...
    'sessionNames', {{}}, ...
    'sessionDates', [], ...
    'nSessions', 0, ...
    'nSessionsTotal', 0, ...
    'nTrialsBySession', [], ...
    'nTrials', 0, ...
    'nTrialsByStep', [], ...
    'scaleBySession', [], ...
    'scaleSEMBySession', [], ...
    'scaleCILoBySession', [], ...
    'scaleCIHiBySession', [], ...
    'prefEnergyBySession', [], ...
    'fitR2BySession', [], ...
    'includeMask', [], ...
    'excludeReasons', {{}}, ...
    'sourceSummaryFiles', {{}}, ...
    'sessionStructs', [] );

end

% ========================================================================
function empirical = computeEmpiricalSummaries(offsetData, opts)

empirical = repmat(struct( ...
    'probeOffsetDeg', NaN, ...
    'stepType', NaN, ...
    'pooledScale', NaN, ...
    'meanSessionScale', NaN, ...
    'medianSessionScale', NaN, ...
    'semSessionScale', NaN, ...
    'boot68', [NaN NaN], ...
    'boot95', [NaN NaN], ...
    'nSessions', 0, ...
    'nTrials', 0, ...
    'nTrialsByStep', [NaN NaN], ...
    'sessionWeighting', 'inverse_variance_kernel_pool', ...
    'distributionSkew', NaN), numel(offsetData), 1);

for k = 1:numel(offsetData)
    x = offsetData(k).scaleBySession;
    empirical(k).probeOffsetDeg   = offsetData(k).probeOffsetDeg;
    empirical(k).stepType         = 1 + strcmp(opts.ScaleStepType, 'inc');
    empirical(k).pooledScale      = computeOffsetPooledScale(offsetData(k), opts, false);
    empirical(k).meanSessionScale = mean(x, 'omitnan');
    empirical(k).medianSessionScale = median(x, 'omitnan');
    empirical(k).semSessionScale  = localSEM(x);
    empirical(k).sessionScale68 = percentileCI(x, 68);
    empirical(k).sessionScale95 = percentileCI(x, 95);
    empirical(k).nSessions        = offsetData(k).nSessions;
    empirical(k).nTrials     = offsetData(k).nTrials;
    empirical(k).nTrialsByStep  = offsetData(k).nTrialsByStep;
    empirical(k).distributionSkew = localSkewness(x);

    % Backward compatibility with earlier field names.
    empirical(k).meanScale   = empirical(k).pooledScale;
    empirical(k).medianScale = empirical(k).medianSessionScale;
    empirical(k).semScale    = empirical(k).semSessionScale;
end

for k = 1:numel(empirical)
    for c = 1:numel(opts.CILevels)
        lvl = opts.CILevels(c);
        fname = sprintf('boot%d', round(lvl));
        empirical(k).(fname) = percentileCI(offsetData(k).scaleBySession, lvl);
    end
    fprintf('      offset %.0f: included %2d/%2d sessions, %5d trials\n', offsetData(k).probeOffsetDeg, ...
                  offsetData(k).nSessions, offsetData(k).nSessionsTotal, offsetData(k).nTrials);
end

end

% ========================================================================
function bootstrap = runHierarchicalBootstrap(offsetData, opts)
% Hierarchical bootstrap over sessions and trials.
%
% Purpose:
%   Estimate uncertainty in the pooled psychophysical scale at each tested
%   probe offset. These offset variances are then used as weights for
%   MT-readout DOG fits.
%
% This function intentionally does not fit any offset-space descriptive
% model. The active model fit is performed later by fitReadoutDOGToScales.

nOffsets = numel(offsetData);
nBoot    = opts.NBoot;

bootScaleMat = nan(nBoot, nOffsets);

% Collect bootstrap pooled scales for each offset.
for b = 1:nBoot
  for k = 1:nOffsets
    bootScaleMat(b, k) = computeOffsetPooledScale(offsetData(k), opts, true);
  end
  if b == 1 || mod(b, max(1, floor(nBoot/10))) == 0
    fprintf('      bootstrap %d / %d\n', b, nBoot);
  end
end

% Estimate offset-specific variances from the bootstrap distribution.
offsetFitVar = nan(1, nOffsets);
for k = 1:nOffsets
  x = bootScaleMat(:, k);
  x = x(isfinite(x));

  if numel(x) >= 2
    offsetFitVar(k) = var(x, 0);
  else
    offsetFitVar(k) = NaN;
  end
end

% Variance floor for stable inverse-variance weighting.
finiteVars = offsetFitVar(isfinite(offsetFitVar) & offsetFitVar > 0);
if isempty(finiteVars)
    varFloor = 1e-6;
else
    varFloor = max(1e-6, 0.01 * median(finiteVars));
end

% Weights for the measured offsets only.
% The 0-deg anchor is not part of offsetData and is imposed analytically by
% the normalized readout prediction S(0)=1.
offsetFitWeights = nan(1, nOffsets);
for k = 1:nOffsets
  if isfinite(offsetFitVar(k)) && offsetFitVar(k) > 0
    offsetFitWeights(k) = 1 ./ max(offsetFitVar(k), varFloor);
  else
    offsetFitWeights(k) = 1;
  end
end

bootstrap = struct();
bootstrap.bootScaleMat = bootScaleMat;

bootstrap.offsetBootstrap = repmat(struct( ...
  'probeOffsetDeg', NaN, ...
  'bootScale', [], ...
  'bootMean', NaN, ...
  'bootMedian', NaN, ...
  'boot68', [NaN NaN], ...
  'boot95', [NaN NaN]), nOffsets, 1);

for k = 1:nOffsets
  x = bootScaleMat(:, k);
  bootstrap.offsetBootstrap(k).probeOffsetDeg = offsetData(k).probeOffsetDeg;
  bootstrap.offsetBootstrap(k).bootScale  = x;
  bootstrap.offsetBootstrap(k).bootMean   = mean(x, 'omitnan');
  bootstrap.offsetBootstrap(k).bootMedian = median(x, 'omitnan');
  bootstrap.offsetBootstrap(k).boot68     = percentileCI(x, 68);
  bootstrap.offsetBootstrap(k).boot95     = percentileCI(x, 95);

  for c = 1:numel(opts.CILevels)
    lvl = opts.CILevels(c);
    bootstrap.offsetBootstrap(k).(sprintf('boot%d', round(lvl))) = ...
      percentileCI(x, lvl);
  end
end

% Keep the name fitBootstrap for minimal downstream disruption, but this is
% now only the bootstrap-derived weighting information for the readout fit.
bootstrap.fitBootstrap = struct( ...
    'modelName', 'none_offset_variance_only', ...
    'offsetsDeg', [offsetData.probeOffsetDeg], ...
    'offsetFitVar', offsetFitVar, ...
    'varFloor', varFloor, ...
    'offsetFitWeights', offsetFitWeights, ...
    'fitWeights', offsetFitWeights, ...
    'note', ['No bootstrap model refits are performed. These fields provide ' ...
             'bootstrap variance estimates and inverse-variance weights for ' ...
             'readout DOG fits.'] );
end

% ========================================================================
function pooledScale = computeOffsetPooledScale(offsetStruct, opts, doBootstrap)
% Mirror kernelAverage: recompute session kernels, normalize the probe
% stream to the pref-noise amplitude convention, pool with inverse-variance
% weighting, then compute scale from the normalized pooled kernels.
sessions = offsetStruct.sessionStructs;
nSess = numel(sessions);
if nSess == 0
    pooledScale = NaN;
    return;
end
if doBootstrap
    sessIdx = randi(nSess, [1 nSess]);
else
    sessIdx = 1:nSess;
end
sessionKernelsNorm = cell(1, numel(sessIdx));
sessionKVarsNorm   = cell(1, numel(sessIdx));
refFrameRateHz = NaN;
refNFrames = NaN;

for j = 1:numel(sessIdx)
    src = sessions(sessIdx(j));

    prefix = '/Users/Shared/Data/IDReadout';
    if ~startsWith(src.noiseFile, [prefix '/IDR'])
      src.noiseFile = [prefix '/IDR' extractAfter(src.noiseFile, strlength(prefix))];
    end

    sessionData = load(src.noiseFile, 'sessionHeader', 'sessionProbeHeader', 'prefNoiseByPatch', 'probeNoiseByPatch', ...
        'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll', 'lr');
    trialIdx = [];
    if doBootstrap
      nTrials = size(sessionData.prefNoiseByPatch, 3);
      trialIdx = randi(nTrials, [1 nTrials]);
    end
    [kernels, kVars] = computeSessionKernels(sessionData, trialIdx);
    [kernelsNorm, kVarsNorm] = normalizeProbeKernelsToPrefAmplitude( ...
      kernels, kVars, src.sessionHeader, src.sessionProbeHeader);

    sessionKernelsNorm{j} = kernelsNorm;
    sessionKVarsNorm{j}   = kVarsNorm;

    thisFrameRateHz = src.sessionHeader.frameRateHz;

    thisNFrames = size(kernels, 4);
    if j == 1
        refFrameRateHz = thisFrameRateHz;
        refNFrames = thisNFrames;
    else
        if thisFrameRateHz ~= refFrameRateHz || thisNFrames ~= refNFrames
            error('Incompatible sessions within offset %.1f deg.', offsetStruct.probeOffsetDeg);
        end
    end
end
[avgKernelsNorm, ~] = poolSessionKernels(sessionKernelsNorm, sessionKVarsNorm, refNFrames);
msPerVFrame = 1000.0 / refFrameRateHz;
[scaleMat, ~, ~, ~] = kernelScaleFit(avgKernelsNorm, msPerVFrame);
pooledScale = selectCompStatsEntry(scaleMat, opts.ScaleSideType, opts.ScaleStepType);
end

% ========================================================================
function history = updateHistory(acrossOffsetSummary, opts)

history = struct( ...
    'runDate', datetime('now'), ...
    'nOffsets', numel(acrossOffsetSummary.empirical), ...
    'offsetsDeg', [acrossOffsetSummary.empirical.probeOffsetDeg], ...
    'nSessionsByOffset', [acrossOffsetSummary.empirical.nSessions], ...
    'primaryFitObject', 'DOG readout over MT preferred direction', ...
    'nFreeParams', acrossOffsetSummary.readoutModel.nFreeParams, ...
    'medianParams', nan, ...
    'saveFile', opts.SaveFile );
end

% ========================================================================
function saveAcrossOffsetSummary(opts, acrossOffsetSummary)

[saveDir, ~, ~] = fileparts(opts.SaveFile);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
save(opts.SaveFile, 'acrossOffsetSummary', '-v7.3');

end

% ========================================================================
function v = selectCompStatsEntry(arr, sideType, stepType)
% Helper for common [sideType x stepType] storage.

sideIdx = mapSideType(sideType);
stepIdx = mapStepType(stepType);

if ndims(arr) < 2
    error('compStats scale array has unexpected dimensions.');
end
v = arr(sideIdx, stepIdx);

end

% ========================================================================
function idx = mapSideType(sideType)

if isnumeric(sideType)
  idx = sideType;
  return;
end

switch lower(char(string(sideType)))
  case 'diff'
    idx = 1;
  case 'change'
    idx = 2;
  case 'nochange'
    idx = 3;
  case 'rf'
    idx = 4;
  case 'opp'
    idx = 5;
  case 'chosen'
    idx = 6;
  case 'notchosen'
    idx = 7;
  otherwise
    error('Unknown side type: %s', char(string(sideType)));
end

end

% ========================================================================
function idx = mapStepType(stepType)

if isnumeric(stepType)
    idx = stepType;
    return;
end

switch lower(char(string(stepType)))
    case 'dec'
        idx = 1;
    case 'inc'
        idx = 2;
    otherwise
        error('Unknown step type: %s', char(string(stepType)));
end

end

% ========================================================================
function se = localSEM(x)

x = x(isfinite(x));
if numel(x) < 2
    se = NaN;
else
    se = std(x) / sqrt(numel(x));
end

end

% ========================================================================
function sk = localSkewness(x)

x = x(isfinite(x));
if numel(x) < 3
    sk = NaN;
else
    sk = skewness(x);
end

end

% ========================================================================
function ci = percentileCI(x, level)

x = x(isfinite(x));
if isempty(x)
    ci = [NaN NaN];
    return;
end
alpha = (100 - level) / 2;
ci = prctile(x, [alpha, 100 - alpha]);

end

% ========================================================================
function [avgKernels, avgKVars] = poolSessionKernels(sessionKernels, sessionKVars, nFrames)
% Pool session kernels using inverse-variance weighting.
%
% Mirrors the helper used by kernelAverage.m.

nSessions = numel(sessionKernels);
nSideTypes = size(sessionKernels{1}, 1);

avgKernels = zeros(nSideTypes, 2, 2, nFrames);
sumWeights = zeros(nSideTypes, 2, 2);

for iSession = 1:nSessions
    kernels = sessionKernels{iSession};
    kVars   = sessionKVars{iSession};

    for sideType = 1:nSideTypes
        for stepType = 1:2
            for streamType = 1:2
                if isfinite(kVars(sideType, stepType, streamType)) && ...
                        kVars(sideType, stepType, streamType) > 0
                    w = 1.0 / kVars(sideType, stepType, streamType);
                    avgKernels(sideType, stepType, streamType, :) = ...
                        avgKernels(sideType, stepType, streamType, :) + ...
                        kernels(sideType, stepType, streamType, :) * w;
                    sumWeights(sideType, stepType, streamType) = ...
                        sumWeights(sideType, stepType, streamType) + w;
                end
            end
        end
    end
end

avgKVars = nan(nSideTypes, 2, 2);
for sideType = 1:nSideTypes
    for stepType = 1:2
        for streamType = 1:2
            if sumWeights(sideType, stepType, streamType) > 0
                avgKernels(sideType, stepType, streamType, :) = ...
                    avgKernels(sideType, stepType, streamType, :) ./ ...
                    sumWeights(sideType, stepType, streamType);
                avgKVars(sideType, stepType, streamType) = ...
                    1.0 / sumWeights(sideType, stepType, streamType);
            else
                avgKernels(sideType, stepType, streamType, :) = nan;
                avgKVars(sideType, stepType, streamType) = nan;
            end
        end
    end
end

end
% ========================================================================
function [kernelsNorm, kVarsNorm, normInfo] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionHeader, sessionProbeHeader)

normInfo = normalizationInfoFromHeaders(sessionHeader, sessionProbeHeader);
probeNormFactor = normInfo.probeNormFactor;

kernelsNorm = kernels;
kVarsNorm   = kVars;

kernelsNorm(:, :, 2, :) = kernelsNorm(:, :, 2, :) * probeNormFactor;
kVarsNorm(:, :, 2)      = kVarsNorm(:, :, 2) * probeNormFactor^2;
end

% ========================================================================
function normInfo = normalizationInfoFromHeaders(sessionHeader, sessionProbeHeader)

prefCohNoisePC  = sessionHeader.prefCohNoisePC;
probeCohNoisePC = sessionProbeHeader.probeCohNoisePC;

normInfo = struct();
normInfo.prefCohNoisePC = prefCohNoisePC;
normInfo.probeCohNoisePC = probeCohNoisePC;

nYokedProbeStreams = probeStreamCountFromSessionProbeHeader(sessionProbeHeader);
combinedProbeCohNoisePC = nYokedProbeStreams * probeCohNoisePC;

normInfo.nYokedProbeStreams = nYokedProbeStreams;
normInfo.combinedProbeCohNoisePC = combinedProbeCohNoisePC;
normInfo.probeNormFactor = (prefCohNoisePC / combinedProbeCohNoisePC)^2;
normInfo.method = ...
  'probe kernels multiplied by (prefCohNoisePC/(nYokedProbeStreams*probeCohNoisePC))^2 before computing normalized ratios/scales';
end

% ========================================================================
function n = probeStreamCountFromSessionProbeHeader(sessionProbeHeader)

assert(isfield(sessionProbeHeader, 'probeDirDeg'), 'updateAcrossOffsetSummaries:MissingProbeDir', ...
    'Cannot determine probe stream count because sessionProbeHeader.probeDirDeg is missing.');

probeDirDeg = abs(double(sessionProbeHeader.probeDirDeg));

if probeDirDeg > 0 && probeDirDeg < 180
    n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
    n = 1;
else
    error('updateAcrossOffsetSummaries:UnsupportedProbeDir', ...
        'Unsupported probeDirDeg for probe normalization: %g.', probeDirDeg);
end
end
