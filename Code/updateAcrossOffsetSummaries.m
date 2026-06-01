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
  dataDir = fullfile(folderPath(), 'Data');
end
opts = parseInputs(dataDir, varargin{:});

% get info about all kernel files meeting criteria
files = {};
fileInfo = table();
for iDir = 1:numel(opts.kernelsDir)
  thisDir = opts.kernelsDir{iDir};
  selectionArgs = [{'FilePattern', opts.FilePattern}, opts.FileSelectionArgs];
  [theFiles, theFileInfo] = selectAnalysisFiles(thisDir, selectionArgs{:});
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
    warning('No usable session summaries found in %s.', summaryDir);
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
assert(numel(bootstrapVarAll) == numel(offsetsDegAll), 'bootstrap offset variances do not match empirical offsets');

isAnchor = abs(offsetsDegAll) < 1e-9;

validFitOffset = ~isAnchor & ...
    nSessionsAll > 0 & ...
    isfinite(offsetsDegAll) & ...
    isfinite(pooledScaleAll) & ...
    isfinite(bootstrapVarAll) & ...
    bootstrapVarAll > 0;
% Fit only non-anchor offsets. The 0-deg value is the normalization anchor. 
fitOffsetsDeg = offsetsDegAll(validFitOffset);
fitScales     = pooledScaleAll(validFitOffset);
fitVars       = bootstrapVarAll(validFitOffset);

% Fixed MT forward model used to map the DOG readout onto predicted scale.
mtModel = makeMTReadoutForwardModel('sigmaMTDeg', 37.5, 'phiDeg', -180:1:179);

% Store measured pooled scales by offset for downstream plotting/fits.
acrossOffsetSummary.measurements = struct();
acrossOffsetSummary.measurements.offsetsDeg   = offsetsDegAll(:)';
acrossOffsetSummary.measurements.pooledScale  = pooledScaleAll(:)';
acrossOffsetSummary.measurements.bootstrapVar = bootstrapVarAll(:)';


readoutModels = struct();
readoutModels.signedDOG = buildReadoutDOGModelSummary( ...
    'signedDOG', 'signed', ...
    fitOffsetsDeg, fitScales, fitVars, ...
    offsetsDegAll, pooledScaleAll, bootstrapVarAll, mtModel, opts);
readoutModels.rectifiedDOG = buildReadoutDOGModelSummary( ...
    'rectifiedDOG', 'rectified', ...
    fitOffsetsDeg, fitScales, fitVars, ...
    offsetsDegAll, pooledScaleAll, bootstrapVarAll, mtModel, opts);

acrossOffsetSummary.readoutModels = readoutModels;
acrossOffsetSummary.readoutModelComparison = compareReadoutDOGModels(readoutModels);

% Backward compatibility: keep the historical top-level readoutModel field
% as the standard signed-template DOG fit.
acrossOffsetSummary.readoutModel = readoutModels.signedDOG;
acrossOffsetSummary.history = updateHistory(acrossOffsetSummary, opts);

saveAcrossOffsetSummary(opts, acrossOffsetSummary);
if opts.MakePlots
  makeAcrossOffsetPlots(acrossOffsetSummary, opts);
end
if opts.Verbose
  fprintf('updateAcrossOffsetSummaries: done.\n');
end

end

% ========================================================================
function rm = buildReadoutDOGModelSummary(modelName, templateMode, ...
    fitOffsetsDeg, fitScales, fitVars, offsetsDegAll, pooledScaleAll, ...
    bootstrapVarAll, mtModel, opts)
% Build one DOG readout model summary for a specified MT template mode.

rm = struct();
rm.activeModelName = modelName;
rm.templateMode = templateMode;
rm.sourceScaleSideType = opts.ScaleSideType;
rm.sourceScaleStepType = opts.ScaleStepType;
rm.sourceScaleMode     = opts.ScaleMode;
rm.mtForwardModelParams = struct( ...
    'sigmaMTDeg', mtModel.sigmaMTDeg, ...
    'phiDeg', mtModel.phiDeg);
rm.measurements = struct();
rm.measurements.offsetsDeg   = offsetsDegAll(:)';
rm.measurements.pooledScale  = pooledScaleAll(:)';
rm.measurements.bootstrapVar = bootstrapVarAll(:)';
rm.readoutNormalization = 'a(0)=1';
rm.note = ['Three-parameter DOG readout fit. No additive baseline is fit, ' ...
           'because any constant readout component is not identifiable ' ...
           'against mean-subtracted MT forward templates.'];

if strcmp(templateMode, 'rectified')
    rm.note = ['Three-parameter DOG readout fit using rectified MT increment ' ...
               'templates. For paired probes, each signed component is ' ...
               'rectified relative to the 0% coherence baseline before the ' ...
               'two components are summed.'];
end
if numel(fitOffsetsDeg) >= 1 && ...
        all(isfinite(fitScales)) && ...
        all(isfinite(fitVars)) && ...
        all(fitVars > 0)
    readoutFit = fitReadoutDOGToScales( ...
        fitOffsetsDeg, fitScales, fitVars, mtModel, ...
        'Bounds', opts.Bounds, ...
        'TemplateMode', templateMode);
    rm.fit = readoutFit;
    rm.nFreeParams = readoutFit.nFreeParams;
    rm.phiDeg = mtModel.phiDeg;
    rm.paramNames = readoutFit.paramNames;
    rm.params = readoutFit.params;
    rm.paramStruct = readoutFit.paramStruct;
    rm.plotOffsetsDeg = 0:1:180;

    hasUsableFit = isfield(readoutFit, 'fitUsable') && readoutFit.fitUsable;

    if hasUsableFit
        rm.readoutPhiRaw = readoutFit.readoutPhiRaw;
        rm.readoutPhi = readoutFit.readoutPhi;

        rm.predictedAtMeasuredOffsets = ...
            predictNormalizedScaleFromReadout(readoutFit.params, offsetsDegAll, mtModel, ...
            'TemplateMode', templateMode);

        rm.plotPredictedScale = ...
            predictNormalizedScaleFromReadout(readoutFit.params, rm.plotOffsetsDeg, mtModel, ...
            'TemplateMode', templateMode);
    else
        rm.readoutPhiRaw = nan(size(mtModel.phiDeg));
        rm.readoutPhi = nan(size(mtModel.phiDeg));
        rm.predictedAtMeasuredOffsets = [];
        rm.plotPredictedScale = [];
    end

else
    rm.fit = [];
    rm.note = ['Readout model not fit: insufficient valid non-anchor offsets ' ...
               'or invalid scale/variance inputs for DOG readout fitting.'];
    rm.nFreeParams = NaN;
    rm.phiDeg = mtModel.phiDeg;
    rm.readoutPhiRaw = nan(size(mtModel.phiDeg));
    rm.readoutPhi = nan(size(mtModel.phiDeg));
    rm.paramNames = {};
    rm.params = [];
    rm.paramStruct = struct();
    rm.predictedAtMeasuredOffsets = [];
    rm.plotOffsetsDeg = 0:1:180;
    rm.plotPredictedScale = [];
end
end

% ========================================================================
function comparison = compareReadoutDOGModels(readoutModels)
% Compact comparison of the signed and rectified DOG fits.

comparison = struct();
comparison.modelNames = {'signedDOG', 'rectifiedDOG'};
comparison.loss = [NaN NaN];
comparison.aic = [NaN NaN];
comparison.aicc = [NaN NaN];
comparison.bic = [NaN NaN];
comparison.reducedChiSq = [NaN NaN];
comparison.pValue = [NaN NaN];
comparison.deltaLoss_signedMinusRectified = NaN;
comparison.deltaAIC_signedMinusRectified = NaN;
comparison.deltaAICc_signedMinusRectified = NaN;
comparison.deltaBIC_signedMinusRectified = NaN;
comparison.preferredByLoss = '';
comparison.preferredByAICc = '';

if isfield(readoutModels, 'signedDOG') && isfield(readoutModels.signedDOG, 'fit') && ...
        ~isempty(readoutModels.signedDOG.fit) && isfield(readoutModels.signedDOG.fit, 'goodnessOfFit')
    g = readoutModels.signedDOG.fit.goodnessOfFit;
    comparison.loss(1) = g.weightedLoss;
    comparison.aic(1) = g.aic;
    comparison.aicc(1) = g.aicc;
    comparison.bic(1) = g.bic;
    comparison.reducedChiSq(1) = g.reducedChiSq;
    comparison.pValue(1) = g.pValue;
end

if isfield(readoutModels, 'rectifiedDOG') && isfield(readoutModels.rectifiedDOG, 'fit') && ...
        ~isempty(readoutModels.rectifiedDOG.fit) && isfield(readoutModels.rectifiedDOG.fit, 'goodnessOfFit')
    g = readoutModels.rectifiedDOG.fit.goodnessOfFit;
    comparison.loss(2) = g.weightedLoss;
    comparison.aic(2) = g.aic;
    comparison.aicc(2) = g.aicc;
    comparison.bic(2) = g.bic;
    comparison.reducedChiSq(2) = g.reducedChiSq;
    comparison.pValue(2) = g.pValue;
end

comparison.deltaLoss_signedMinusRectified = comparison.loss(1) - comparison.loss(2);
comparison.deltaAIC_signedMinusRectified  = comparison.aic(1)  - comparison.aic(2);
comparison.deltaAICc_signedMinusRectified = comparison.aicc(1) - comparison.aicc(2);
comparison.deltaBIC_signedMinusRectified  = comparison.bic(1)  - comparison.bic(2);

if all(isfinite(comparison.loss))
    if comparison.loss(1) <= comparison.loss(2)
        comparison.preferredByLoss = 'signedDOG';
    else
        comparison.preferredByLoss = 'rectifiedDOG';
    end
end
if all(isfinite(comparison.aicc))
    if comparison.aicc(1) <= comparison.aicc(2)
        comparison.preferredByAICc = 'signedDOG';
    else
        comparison.preferredByAICc = 'rectifiedDOG';
    end
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

addRequired(p, 'dataDir', @(x) ischar(x) || isstring(x) || iscell(x));
defaultSaveFile = fullfile(folderPath(), 'Data', 'AcrossOffsetSummaries', 'IDR_acrossOffsetSummary.mat');
addParameter(p, 'SaveFile', defaultSaveFile, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotDir',  fullfile(dataFolderPath(), '..', 'Plots', 'ReadoutFits'), @(x) ischar(x) || isstring(x));
addParameter(p, 'NBoot', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CILevels', [68 95], @(x) isnumeric(x) && isvector(x) && all(x > 0) && all(x < 100));
addParameter(p, 'Bin179With180', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'FileSelectionArgs', {}, @(x) iscell(x));
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
opts.SaveFile = char(opts.SaveFile);
opts.PlotDir  = char(opts.PlotDir);
opts.FilePattern = char(opts.FilePattern);
opts.OffsetField = char(opts.OffsetField);
opts.SessionNameField = char(opts.SessionNameField);
opts.ScaleMode = char(opts.ScaleMode);


if ~exist(fileparts(opts.SaveFile), 'dir')
  mkdir(fileparts(opts.SaveFile));
end
if opts.MakePlots && ~exist(opts.PlotDir, 'dir')
  mkdir(opts.PlotDir);
end

%inset a convenience flag for whether we are binning 179° with 180°
names = opts.FileSelectionArgs(1:2:end);
values = opts.FileSelectionArgs(2:2:end);
isMatch = cellfun(@(x) (ischar(x) || isstring(x)) && strcmpi(char(x), 'Bin179With180'), names);
if any(isMatch)
  opts.Bin179With180 = values{find(isMatch, 1, 'last')};
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

    S = load(thisFile);
    sessionRecord = emptySessionStruct();
    sessionRecord.sessionHeader = S.sessionHeader;
    sessionRecord.sessionProbeHeader = S.sessionProbeHeader;
    sessionRecord.probeOffsetDeg = S.sessionProbeHeader.probeDirDeg;
    sessionRecord.sessionName = S.sessionHeader.fileName;
    sessionRecord.kernelFile = S.sessionProbeHeader.kernelFile;
    sessionRecord.noiseFile  = S.sessionProbeHeader.noiseFile;
    sessionRecord.compStats = S.compStats;
    sessionRecord.hitStats  = S.hitStats;
    sessionRecord.sourceFile = thisFile;
    sessionRecord.scalePointEstimate = selectCompStatsEntry(S.compStats.normScale, opts.ScaleSideType, opts.ScaleStepType);
    sessionRecord.nTrials = sum(S.hitStats.nTrials);
    sessionRecord.nTrialsByStep = S.hitStats.nTrials;
    if ~isempty(fileInfo) && height(fileInfo) >= iFile
        sessionRecord.fileInfo = fileInfo(iFile, :);
    end

    [tfExclude, reasonStr] = excludeFromReadoutFit(sessionRecord.probeOffsetDeg, ...
              sessionRecord.sessionProbeHeader.probeCohNoisePC, opts.Bin179With180);
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
    'kernelFile', '', ...
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

function [exclude, reason] = excludeFromReadoutFit(probeDirDeg, probeCohNoisePC, Bin179With180)
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
if Bin179With180
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

% if we are binning 179° with 180°, combine them if both are preset
if opts.Bin179With180
  index179 = find(offsetKeys == 179, 1);
  index180 = find(offsetKeys == 180, 1);
  if ~isempty(index179) && ~isempty(index180)
    offsetData(index179).sessionNames = [offsetData(index179).sessionNames, offsetData(index180).sessionNames];
    offsetData(index179).sessionDates = [offsetData(index179).sessionDates, offsetData(index180).sessionDates];
    offsetData(index179).nSessions = offsetData(index179).nSessions + offsetData(index180).nSessions;
    offsetData(index179).nSessionsTotal = offsetData(index179).nSessionsTotal + offsetData(index180).nSessionsTotal;
    offsetData(index179).nTrialsBySession = [offsetData(index179).nTrialsBySession, offsetData(index180).nTrialsBySession];
    offsetData(index179).nTrialsByStep = ...
          offsetData(index179).nTrialsByStep + offsetData(index180).nTrialsByStep;
    offsetData(index179).nTrials = offsetData(index179).nSessionsTotal + offsetData(index180).nTrials;
    offsetData(index179).scaleBySession = [offsetData(index179).scaleBySession, offsetData(index180).scaleBySession];

    offsetData(index179).scaleSEMBySession = nan(size(offsetData(index179).scaleBySession));
    offsetData(index179).scaleSEMBySession = nan(size(offsetData(index179).scaleBySession));
    offsetData(index179).scaleCILoBySession = nan(size(offsetData(index179).scaleBySession));
    offsetData(index179).scaleCIHiBySession = nan(size(offsetData(index179).scaleBySession));
    offsetData(index179).prefEnergyBySession = nan(size(offsetData(index179).scaleBySession));
    offsetData(index179).fitR2BySession = nan(size(offsetData(index179).scaleBySession));

    offsetData(index179).includeMask = [offsetData(index179).includeMask, offsetData(index180).includeMask];
    offsetData(index179).excludeReasons = [offsetData(index179).excludeReasons, offsetData(index180).excludeReasons];
    offsetData(index179).sourceSummaryFiles = [offsetData(index179).sourceSummaryFiles, offsetData(index180).sourceSummaryFiles];
    offsetData(index179).sessionStructs = [offsetData(index179).sessionStructs; offsetData(index180).sessionStructs];

    offsetData(index180) = [];
    offsetKeys(index180) = [];
  end
end

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
    fprintf('      offset %.1f: included %d/%d sessions, %d trials\n', offsetData(k).probeOffsetDeg, ...
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
function scaleVal = recomputeSessionScalePointEstimate(sessionRecord, opts)
% Recompute the requested side/step point estimate from source noise matrices.

if isempty(sessionRecord.noiseFile) || ~exist(sessionRecord.noiseFile, 'file')
    scaleVal = sessionRecord.scalePointEstimate;
    return;
end

sessionData = load(sessionRecord.noiseFile);
[~, ~, ~, ~, compStats] = computeSessionKernels(sessionData, []);

if isfield(compStats, 'normScale')
    scaleVal = selectCompStatsEntry(compStats.normScale, opts.ScaleSideType, opts.ScaleStepType);
else
    warning('updateAcrossOffsetSummaries:MissingNormScale', ...
        'compStats.normScale missing for session %s; falling back to raw scale.', ...
        sessionRecord.sessionName);
    scaleVal = selectCompStatsEntry(compStats.scale, opts.ScaleSideType, opts.ScaleStepType);
end
end

% % ========================================================================
% function [kernels, kVars, hitStats, compStats] = recomputeSessionKernelStruct(sessionRecord, ~, doTrialBootstrap)
% % Recompute full session kernel outputs from source noise matrices.
% 
% sessionData = load(sessionRecord.noiseFile);
% trialIdx = [];
% if doTrialBootstrap
%     nTrials = size(sessionData.prefNoiseByPatch, 3);
%     trialIdx = randi(nTrials, [1 nTrials]);
% end
% [kernels, kVars, ~, hitStats, compStats] = computeSessionKernels(sessionData, trialIdx);
% end

% ========================================================================
function deltaM = computeMTSymmetrizedDeltaM(phiDeg, deltaDeg, mtModel)
% Compute the effective MT template for the probe-noise configuration.
%
% Conventions:
%   - delta = 0 deg: single noise channel
%   - 0 < delta < 180 deg: symmetric probe noise at +delta and -delta.
%     Each stream has amplitude 1/sqrt(2) relative to the single-stream
%     probe amplitude, so the effective template is
%         (1/sqrt(2)) * [Delta m(+delta) + Delta m(-delta)].
%   - delta = 180 deg: single noise channel

phiDeg = phiDeg(:)';
if numel(phiDeg) ~= numel(mtModel.phiDeg) || any(phiDeg ~= mtModel.phiDeg)
    error('computeMTSymmetrizedDeltaM requires the shared phi-grid from mtModel.phiDeg.');
end

deltaDeg = abs(deltaDeg);

if abs(deltaDeg) < 1e-9 || abs(deltaDeg - 180) < 1e-9
    deltaM = computeMTDeltaM(phiDeg, deltaDeg, mtModel);
else
    deltaMPlus  = computeMTDeltaM(phiDeg, +deltaDeg, mtModel);
    deltaMMinus = computeMTDeltaM(phiDeg, -deltaDeg, mtModel);
    deltaM = 0.5 .* (deltaMPlus + deltaMMinus);
end
end

% ========================================================================
function deltaM = computeMTSymmetrizedTemplate(phiDeg, deltaDeg, mtModel, templateMode)
% Return the effective MT template for the requested forward-model mode.
%
% signed:
%   Use the historical signed, mean-subtracted MT perturbation.
% rectified:
%   Half-wave rectify each single-stream signed perturbation relative to the
%   0% coherence baseline before summing paired +delta/-delta components.

if nargin < 4 || isempty(templateMode)
    templateMode = 'signed';
end
templateMode = lower(char(string(templateMode)));

phiDeg = phiDeg(:)';
if numel(phiDeg) ~= numel(mtModel.phiDeg) || any(phiDeg ~= mtModel.phiDeg)
    error('computeMTSymmetrizedTemplate requires the shared phi-grid from mtModel.phiDeg.');
end

deltaDeg = abs(deltaDeg);

switch templateMode
    case 'signed'
        deltaM = computeMTSymmetrizedDeltaM(phiDeg, deltaDeg, mtModel);

    case 'rectified'
        if abs(deltaDeg) < 1e-9 || abs(deltaDeg - 180) < 1e-9
            deltaM = max(computeMTDeltaM(phiDeg, deltaDeg, mtModel), 0);
        else
            deltaMPlus  = max(computeMTDeltaM(phiDeg, +deltaDeg, mtModel), 0);
            deltaMMinus = max(computeMTDeltaM(phiDeg, -deltaDeg, mtModel), 0);
            deltaM = 0.5 .* (deltaMPlus + deltaMMinus);
        end

    otherwise
        error('Unknown MT template mode: %s', templateMode);
end
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
    sessionData = load(src.noiseFile);
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
function mtModel = makeMTReadoutForwardModel(varargin)
% Fixed MT forward model for mapping readout a(phi) onto predicted scale.
%
% Conventions:
%   - 1 deg phi grid
%   - do not double-count -180 and 180
%   - same grid used for readout, MT tuning, and discrete mean subtraction

p = inputParser;
addParameter(p, 'sigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'phiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
parse(p, varargin{:});

phiDeg = p.Results.phiDeg(:)';
sigmaMTDeg = p.Results.sigmaMTDeg;

mtModel = struct();
mtModel.sigmaMTDeg = sigmaMTDeg;
mtModel.phiDeg = phiDeg;
mtModel.gridStepDeg = median(diff(phiDeg));

% Canonical MT tuning template centered at 0 deg.
G0 = exp(-(phiDeg.^2) ./ (2 * sigmaMTDeg.^2));
Gbar = mean(G0);

mtModel.G0 = G0;
mtModel.Gbar = Gbar;
mtModel.note = ['Delta m(phi;delta) is defined as G(phi-delta) minus the ' ...
    'discrete mean of G on the same phi-grid.'];
end

% ========================================================================
function [aPhi, paramStruct] = evaluateReadoutDOG(phiDeg, params)
% Evaluate DOG readout over MT preferred direction.
%
% Three-parameter form, used for the signed-template fit:
%   a(phi) = exp(-(phi^2)/(2*sigmaC^2)) ...
%          - As * exp(-(phi^2)/(2*sigmaS^2))
%
% Four-parameter form, used for the rectified-template fit:
%   a(phi) = B ...
%          + exp(-(phi^2)/(2*sigmaC^2)) ...
%          - As * exp(-(phi^2)/(2*sigmaS^2))
%
% params:
%   [sigmaC, sigmaS, As]       signed-template DOG
%   [sigmaC, sigmaS, As, B]    rectified-template DOG with nonzero asymptote

phiDeg = phiDeg(:)';

if numel(params) == 3
    baselineOffset = 0;
elseif numel(params) == 4
    baselineOffset = params(4);
else
    error(['evaluateReadoutDOG requires params = [sigmaC, sigmaS, As] ' ...
           'or [sigmaC, sigmaS, As, baselineOffset].']);
end

sigmaC = params(1);
sigmaS = params(2);
As     = params(3);

aPhi = baselineOffset ...
     + exp(-(phiDeg.^2) ./ (2 * sigmaC.^2)) ...
     - As .* exp(-(phiDeg.^2) ./ (2 * sigmaS.^2));

paramStruct = struct( ...
    'sigmaCenterDeg', sigmaC, ...
    'sigmaSurroundDeg', sigmaS, ...
    'surroundGain', As, ...
    'baselineOffset', baselineOffset );
end

% ========================================================================
function aPhiNorm = normalizeReadoutAtPreferred(phiDeg, aPhi)
% Normalize readout so that the preferred-direction value a(0) == 1.
%
% This is a convention for identifiability and interpretation. It does not
% affect the predicted normalized scale because the scale prediction is
% invariant to multiplicative rescaling of the readout.

phiDeg = phiDeg(:)';
aPhi   = aPhi(:)';

idx0 = find(abs(phiDeg) < 1e-9, 1, 'first');
if isempty(idx0)
    error('normalizeReadoutAtPreferred requires phiDeg to contain 0 deg.');
end

a0 = aPhi(idx0);
if ~isfinite(a0) || abs(a0) < 1e-12
    error('Cannot normalize readout at preferred direction because a(0) is zero or non-finite.');
end

aPhiNorm = aPhi ./ a0;
end

% ========================================================================
function deltaM = computeMTDeltaM(phiDeg, deltaDeg, mtModel)
% Compute Delta m(phi;delta) on the shared phi-grid.

phiDeg = phiDeg(:)';
if numel(phiDeg) ~= numel(mtModel.phiDeg) || any(phiDeg ~= mtModel.phiDeg)
    error('computeMTDeltaM requires the shared phi-grid from mtModel.phiDeg.');
end

shifted = wrapTo180Local(phiDeg - deltaDeg);
G = exp(-(shifted.^2) ./ (2 * mtModel.sigmaMTDeg.^2));
deltaM = G - mean(G);
end

% ========================================================================
function predScale = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel, varargin)
% Predict normalized scale from a DOG readout and selected MT template mode.
%
% TemplateMode = 'signed' reproduces the historical model:
%   Delta m_eff(phi;delta) = (1/sqrt(2)) * [Delta m(phi;+delta) + Delta m(phi;-delta)]
%
% TemplateMode = 'rectified' rectifies each single-stream signed perturbation
% relative to the 0% coherence baseline before summing paired probes:
%   Delta m_eff(phi;delta) = (1/sqrt(2)) * [max(Delta m(phi;+delta),0) + max(Delta m(phi;-delta),0)]
%
% S_pred(delta) = sum_phi a(phi) Delta m_eff(phi;delta) / ...
%                 sum_phi a(phi) Delta m_eff(phi;0)

p = inputParser;
addParameter(p, 'TemplateMode', 'signed', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
templateMode = lower(char(string(p.Results.TemplateMode)));

phiDeg = mtModel.phiDeg;
[aPhi, ~] = evaluateReadoutDOG(phiDeg, params);

refVal = sum(aPhi .* computeMTSymmetrizedTemplate(phiDeg, 0, mtModel, templateMode));

offsetsDeg = offsetsDeg(:)';
predScale = nan(size(offsetsDeg));
for i = 1:numel(offsetsDeg)
    probeVal = sum(aPhi .* computeMTSymmetrizedTemplate(phiDeg, offsetsDeg(i), mtModel, templateMode));
    if isfinite(refVal) && abs(refVal) > 0
        predScale(i) = probeVal / refVal;
    else
        predScale(i) = NaN;
    end
end
end

% ========================================================================
function plotReadoutDiagnostics(figNum, acrossOffsetSummary, opts)
% Plot fitted readout, MT templates, and their products to visualize how
% overlap determines predicted normalized scale.

  rm = acrossOffsetSummary.readoutModel;
  if isfield(rm, 'templateMode')
      templateMode = rm.templateMode;
  else
      templateMode = 'signed';
  end
  if ~isfield(rm, 'fit') || isempty(rm.fit) || ~rm.fit.fitSuccess
      return;
  end

  phiDeg = rm.phiDeg(:)';
  aPhi   = rm.readoutPhi(:)';   % normalized display readout, a(0)=1
  mtp = rm.mtForwardModelParams;
  mtModel = makeMTReadoutForwardModel('sigmaMTDeg', mtp.sigmaMTDeg, 'phiDeg', mtp.phiDeg);
  offsetsDeg = [0, rm.fit.offsetsDeg];  
  nOffsets = numel(offsetsDeg);
  
  deltaM   = cell(1, nOffsets);
  prodTerm = cell(1, nOffsets);
  overlap  = nan(1, nOffsets);
  posPart  = nan(1, nOffsets);
  negPart  = nan(1, nOffsets);

  for i = 1:nOffsets
      deltaM{i} = computeMTSymmetrizedTemplate(phiDeg, offsetsDeg(i), mtModel, templateMode);
      prodTerm{i} = aPhi .* deltaM{i};
      overlap(i) = sum(prodTerm{i});
      posPart(i) = sum(max(prodTerm{i}, 0));
      negPart(i) = sum(min(prodTerm{i}, 0));
  end
  
  idx0 = find(abs(offsetsDeg) < 1e-9, 1, 'first');
  if isempty(idx0)
      warning('plotReadoutDiagnostics: OffsetsDeg does not include 0. Ratios will not be shown.');
  end
  
  % ---- Build fit-vs-flat comparison figure ----
  prodTermFit  = cell(1, nOffsets);
  prodTermFlat = cell(1, nOffsets);
  for i = 1:nOffsets
      prodTermFit{i}  = aPhi .* deltaM{i};
      prodTermFlat{i} = ones(size(aPhi)) .* deltaM{i};   % flat readout = 1
  end

  % -- set up figure to plot three panels
  fig = figure(figNum); clf;
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
  lineCol = lines(nOffsets);
  
  % ---- top panel: readout function ----
  nexttile; hold on;
  hFitReadout  = plot(phiDeg, aPhi, 'k-', 'LineWidth', 2);
  % hFlatReadout = plot(phiDeg, ones(size(phiDeg)), 'k--', 'LineWidth', 1.5);
  plot(phiDeg, zeros(size(phiDeg)), 'k:', 'LineWidth', 1);
  xlabel('\phi (deg)');
  ylabel('a(\phi)');
  title(sprintf('Fitted DOG readout (%s template, %d bootstraps)', templateMode, opts.NBoot));
  legend(hFitReadout, {'Fitted readout a(\phi)'}, 'Location', 'northeast');
  paramText = cell(rm.nFreeParams, 1);
  for p = 1:rm.nFreeParams
    paramText{p} = sprintf('%s: %.4f', rm.paramNames{p}, rm.params(p));
  end
  text(-100, 0.95, paramText, 'horizontalAlignment', 'right', 'VerticalAlignment', 'top');
  ylimits = ylim();
  ylim([min(0.2, ylimits(1)), max(1.1, ylimits(2))]);
  box off;
  
  % ---- middle panel: MT populations responses ----
  nexttile; hold on;
  title(sprintf('MT Population Responses to Probes (%s template; Flat Readout)', templateMode));
  % MT templates, same colors used in the lower panel
  hTemplates = gobjects(1, nOffsets);
  for i = 1:nOffsets
      hTemplates(i) = plot(phiDeg, deltaM{i}, '-', 'Color', lineCol(i,:), 'LineWidth', 1.5);
  end
  yline(0, 'k:');
  xlabel('\phi (deg)');
  ylabel('\Delta m(\phi;\delta)');
  legend(hTemplates, [arrayfun(@(d) sprintf('\\Delta m(\\phi;%g^\\circ)', d), offsetsDeg, ...
       'UniformOutput', false)], 'Location', 'best');
  box off;
  
  % ---- Bottom panel: overlap contribution functions ----
  nexttile; hold on;
  title('Weighted Population Responses (Fit)');
  legendHandles = gobjects(0);
  legendLabels = {};
  
  for i = 1:nOffsets
      h1 = plot(phiDeg, prodTermFit{i}, '-', 'Color', lineCol(i,:), 'LineWidth', 2);
      legendHandles(end+1) = h1; %#ok<AGROW>
      legendLabels{end+1} = sprintf('%g^\\circ fit, <a,\\Deltam> = %.2g; S_{fit} %.2f', ...
        offsetsDeg(i), overlap(i), overlap(i)/overlap(1)); %#ok<AGROW>
  end
  yline(0, 'k:');
  xlabel('\phi (deg)');
  ylabel('a(\phi)\Delta m(\phi;\delta)');
  legend(legendHandles, legendLabels, 'Location', 'best');
  box off;
  if ~isempty(char(string(opts.SaveFile)))
      saveas(fig, fullfile(opts.PlotDir, sprintf('ReadoutFunctions_%s.pdf', templateMode)));
  end
end

% ========================================================================
function ang = wrapTo180Local(ang)
% Map angles to [-180, 180).

ang = mod(ang + 180, 360) - 180;
end

% ========================================================================
function fitResult = fitReadoutDOGToScales(offsetsDeg, obsScale, obsVar, mtModel, varargin)
% Fit DOG readout to observed non-anchor scale values.
%
% Supports:
%   signed template:    [sigmaCenterDeg, sigmaSurroundDeg, surroundGain]
%   rectified template: [sigmaCenterDeg, sigmaSurroundDeg, surroundGain, baselineOffset]
%
% Weighted objective:
%   sum_i (obsScale_i - predScale_i)^2 / obsVar_i
%
% Uses multi-start fmincon to reduce sensitivity to local minima.

p = inputParser;
addParameter(p, 'Bounds', struct(), @(x) isstruct(x));
addParameter(p, 'TemplateMode', 'signed', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opts = p.Results;
templateMode = lower(char(string(opts.TemplateMode)));

[p0, lb, ub, paramNames] = initialGuessReadoutDOG(offsetsDeg, opts.Bounds, templateMode);

p0 = p0(:).';
lb = lb(:).';
ub = ub(:).';
p0 = min(max(p0, lb), ub);

obj = @(params) readoutDOGObjective(params, offsetsDeg, obsScale, obsVar, mtModel, templateMode);

% ---- Build multi-start list ----
% First row is the historical single-start initial guess. Subsequent rows
% deliberately sample narrow, medium, and broad DOGs. The baseline column is
% appended only for the rectified 4-parameter fit.
baseStarts3 = [
    p0(1:min(3,numel(p0)))
    % 0.01   0.02  2.0
   %  1.5    2.5    0.5
   %  2.5    4.0    0.5
   %  4.0    6.0    0.5
   %  6.0    9.0    0.5
   %  8.0   12.0    0.5
   % 10.0   15.0    0.5
   % 15.0   25.0    0.5
   % 25.0   50.0    0.5
   % 40.0   90.0    0.5
   % 75.0  150.0    0.5
   %  2.5    4.0    1.0
   %  5.0    8.0    1.0
   % 10.0   15.0    1.0
   % 15.0   25.0    1.0
   % 25.0   50.0    1.0
   % 40.0   90.0    1.0
   %  5.0   20.0    0.25
   % 10.0   40.0    0.25
   % 25.0   90.0    0.25
];

% Remove accidental duplicate rows before adding the rectified baseline.
baseStarts3 = unique(baseStarts3, 'rows', 'stable');

if numel(p0) == 3
    startParams = baseStarts3;
elseif numel(p0) == 4
    % Try several baseline offsets. Include zero first so the original
    % rectified behavior remains one of the candidate starts.
    baselineStarts = [0, -0.2, 0.2, -0.5, 0.5];
    startParams = [];
    for iB = 1:numel(baselineStarts)
        startParams = [startParams; ...
            baseStarts3, repmat(baselineStarts(iB), size(baseStarts3, 1), 1)]; %#ok<AGROW>
    end
else
    error('fitReadoutDOGToScales:BadParamCount', ...
        'Expected 3 or 4 DOG parameters, got %d.', numel(p0));
end

% Respect user-supplied bounds and remove starts that cannot be finite.
for iStart = 1:size(startParams, 1)
    startParams(iStart,:) = min(max(startParams(iStart,:), lb), ub);
end
startParams = unique(startParams, 'rows', 'stable');

% ---- Run fmincon from each start and keep the best finite objective ----
bestParams = nan(size(p0));
bestLoss = Inf;
bestExitflag = NaN;
bestStartIndex = NaN;
fitSuccessAny = false;

startLog = repmat(struct( ...
    'startIndex', NaN, ...
    'pStart', [], ...
    'pFit', [], ...
    'loss', NaN, ...
    'exitflag', NaN, ...
    'success', false, ...
    'message', ''), size(startParams, 1), 1);

try
    fminconOpts = optimoptions( ...
        'fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point');

    for iStart = 1:size(startParams, 1)
        pStart = startParams(iStart,:);

        startLog(iStart).startIndex = iStart;
        startLog(iStart).pStart = pStart;

        try
            [pFit, fval, exitflag] = fmincon( ...
                obj, ...
                pStart, ...
                [], [], [], [], ...
                lb, ub, ...
                [], ...
                fminconOpts);

            thisSuccess = exitflag > 0 && all(isfinite(pFit)) && isfinite(fval);

            startLog(iStart).pFit = pFit;
            startLog(iStart).loss = fval;
            startLog(iStart).exitflag = exitflag;
            startLog(iStart).success = thisSuccess;

            if thisSuccess
                fitSuccessAny = true;
            end

            % Keep the best finite solution, even if fmincon reports a
            % nonpositive exit flag. This avoids discarding useful boundary
            % solutions while preserving fitSuccessAny separately.
            if all(isfinite(pFit)) && isfinite(fval) && fval < bestLoss
                bestParams = pFit;
                bestLoss = fval;
                bestExitflag = exitflag;
                bestStartIndex = iStart;
            end

        catch MEstart
            startLog(iStart).message = MEstart.message;
        end
    end

catch ME
    warning('fmincon setup failed in fitReadoutDOGToScales: %s', ME.message);
end

params = bestParams;
loss = bestLoss;
exitflag = bestExitflag;

fitConverged = fitSuccessAny && all(isfinite(params)) && isfinite(loss);
fitUsable = all(isfinite(params)) && isfinite(loss);

predMeasured = nan(size(offsetsDeg(:)'));
paramStruct = struct();

if fitUsable
    predMeasured = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel, ...
        'TemplateMode', templateMode);
    predMeasured = predMeasured(:).';

    fitUsable = all(isfinite(predMeasured));

    if fitUsable
        [~, paramStruct] = evaluateReadoutDOG(mtModel.phiDeg, params);
    end
end

gof = computeReadoutGoodnessOfFit(obsScale, predMeasured, obsVar, numel(p0));
fitUsable = fitUsable && isfinite(gof.weightedLoss);

fitResult = struct();
fitResult.modelName = 'dog_readout';
fitResult.templateMode = templateMode;

fitResult.fitSuccess = fitConverged;
fitResult.loss = loss;
fitResult.exitflag = exitflag;
fitResult.fitConverged = fitConverged;
fitResult.fitUsable = fitUsable;

fitResult.params = params;
fitResult.paramStruct = paramStruct;
fitResult.paramNames = paramNames;
fitResult.nFreeParams = numel(paramNames);

fitResult.offsetsDeg = offsetsDeg(:)';
fitResult.observedScale = obsScale(:)';
fitResult.observedVar = obsVar(:)';
fitResult.predictedScale = predMeasured(:)';
fitResult.residuals = gof.residuals;
fitResult.standardizedResiduals = gof.standardizedResiduals;
fitResult.goodnessOfFit = gof;
fitResult.phiDeg = mtModel.phiDeg;

% New diagnostics for the multi-start search.
fitResult.multistart = struct();
fitResult.multistart.nStarts = size(startParams, 1);
fitResult.multistart.starts = startParams;
fitResult.multistart.bestStartIndex = bestStartIndex;
fitResult.multistart.startLog = startLog;
fitResult.multistart.note = ...
    'Best finite weighted-loss solution retained across multiple fmincon starts.';

if fitUsable
    fitResult.readoutPhiRaw = evaluateReadoutDOG(mtModel.phiDeg, params);
    fitResult.readoutPhi = normalizeReadoutAtPreferred(mtModel.phiDeg, fitResult.readoutPhiRaw);
else
    fitResult.readoutPhiRaw = nan(size(mtModel.phiDeg));
    fitResult.readoutPhi = nan(size(mtModel.phiDeg));
end

end

% ========================================================================
function [p0, lb, ub, paramNames] = initialGuessReadoutDOG(offsetsDeg, bounds, templateMode)

if nargin < 3 || isempty(templateMode)
    templateMode = 'signed';
end
templateMode = char(string(templateMode));

sigmaC0 = max(10, min(50, median(offsetsDeg(offsetsDeg > 0), 'omitnan')));
if isempty(sigmaC0) || ~isfinite(sigmaC0)
    sigmaC0 = 25;
end

sigmaS0 = max(sigmaC0 + 20, 90);
As0 = 0.5;

p0 = [sigmaC0, sigmaS0, As0];
lb = [1e-3, 1e-3, 0];
ub = [300, 300, 10];
paramNames = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};

% For the rectified MT template, the template no longer sums to zero. A
% constant component in the DOG readout is therefore identifiable and lets
% the readout asymptote to a nonzero value at far-from-preferred directions.
% Keep the signed-template fit at three parameters for backward compatibility
% and because this baseline is nulled by the mean-subtracted signed template.
if strcmpi(templateMode, 'rectified')
    baselineOffset0 = 0;
    p0 = [p0, baselineOffset0];
    lb = [lb, -2];
    ub = [ub,  2];
    paramNames = [paramNames, {'baselineOffset'}];
end

if isfield(bounds, 'readoutDOG')
    B = bounds.readoutDOG;
elseif isfield(bounds, 'dog_readout')
    B = bounds.dog_readout;
elseif isfield(bounds, 'dog')
    B = bounds.dog;
else
    B = struct();
end

if isfield(B, 'lb')
    lbUser = B.lb(:)';
    if numel(lbUser) == numel(lb)
        lb = lbUser;
    elseif numel(lbUser) == 3 && numel(lb) == 4
        lb(1:3) = lbUser;
    else
        error('readoutDOG lower bounds must have %d entries, or 3 entries for shared DOG parameters.', numel(lb));
    end
end
if isfield(B, 'ub')
    ubUser = B.ub(:)';
    if numel(ubUser) == numel(ub)
        ub = ubUser;
    elseif numel(ubUser) == 3 && numel(ub) == 4
        ub(1:3) = ubUser;
    else
        error('readoutDOG upper bounds must have %d entries, or 3 entries for shared DOG parameters.', numel(ub));
    end
end

p0 = p0(:);
lb = lb(:);
ub = ub(:);
end

% ========================================================================
function sse = readoutDOGObjective(params, offsetsDeg, obsScale, obsVar, mtModel, templateMode)

predScale = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel, ...
    'TemplateMode', templateMode);
resid = obsScale(:)' - predScale(:)';

valid = isfinite(resid) & isfinite(obsVar(:)') & obsVar(:)' > 0;

if ~any(valid)
    sse = Inf;
    return;
end

sse = sum((resid(valid) .^ 2) ./ obsVar(valid));

% Soft constraint: surround broader than center.
sigmaC = params(1);
sigmaS = params(2);
minRatio = 1.25;
if sigmaS < minRatio * sigmaC
    sse = sse + 1e6 + 1e3 * (minRatio * sigmaC - sigmaS);
end

if ~isfinite(sse)
    sse = Inf;
end
end

% ========================================================================
function gof = computeReadoutGoodnessOfFit(obsScale, predScale, obsVar, nFreeParams)
% Goodness-of-fit diagnostics for weighted DOG scale fits.

obsScale = obsScale(:)';
predScale = predScale(:)';
obsVar = obsVar(:)';

resid = obsScale - predScale;
valid = isfinite(resid) & isfinite(obsVar) & obsVar > 0;

weightedTerms = nan(size(resid));
standardizedResiduals = nan(size(resid));
if any(valid)
    weightedTerms(valid) = (resid(valid).^2) ./ obsVar(valid);
    standardizedResiduals(valid) = resid(valid) ./ sqrt(obsVar(valid));
end

nData = sum(valid);
df = nData - nFreeParams;
if nData > 0
    weightedLoss = sum(weightedTerms(valid));
else
    weightedLoss = NaN;
end

if df > 0 && isfinite(weightedLoss)
    reducedChiSq = weightedLoss / df;
    pValue = gammainc(weightedLoss / 2, df / 2, 'upper');
else
    reducedChiSq = NaN;
    pValue = NaN;
end

% These information criteria use the weighted residual loss without the
% Gaussian normalization constants. That is sufficient for comparing models
% fit to the same observations and variances.
aic = weightedLoss + 2 * nFreeParams;
if nData > nFreeParams + 1
    aicc = aic + (2 * nFreeParams * (nFreeParams + 1)) / (nData - nFreeParams - 1);
else
    aicc = NaN;
end
if nData > 0
    bic = weightedLoss + nFreeParams * log(nData);
else
    bic = NaN;
end

gof = struct();
gof.weightedLoss = weightedLoss;
gof.nData = nData;
gof.nFreeParams = nFreeParams;
gof.df = df;
gof.reducedChiSq = reducedChiSq;
gof.pValue = pValue;
gof.aic = aic;
gof.aicc = aicc;
gof.bic = bic;
gof.residuals = resid;
gof.standardizedResiduals = standardizedResiduals;
gof.weightedLossTerms = weightedTerms;
gof.note = ['pValue is an approximate chi-square tail probability using ' ...
            'bootstrap variances as known observation variances. AIC/AICc/BIC ' ...
            'omit constants common to models fit to the same observations.'];
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
function makeAcrossOffsetPlots(acrossOffsetSummary, opts)
% Primary plots for the DOG readout models:
%   1) observed scale values by probe offset, with signed and rectified DOG
%      predictions overlaid when fits are available
%   2) fitted readout/template diagnostics for each successful model

emp = acrossOffsetSummary.empirical;
offsets = [emp.probeOffsetDeg];
obsScale = [emp.pooledScale];
ci95 = vertcat(acrossOffsetSummary.bootstrap.offsetBootstrap.boot95);

% ---- Plot 1: observed and fit scale by offset ----
fig1 = figure(300); clf; hold on;
hObs = errorbar([0, offsets], [1, obsScale], [0, obsScale - ci95(:,1)'], [0, ci95(:,2)' - obsScale], ...
    'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'k');
plot([0, 180], [0,0], 'k:');
legendHandles = hObs;
legendLabels = {'Observed Scale (95% CI)'};
signedRM = acrossOffsetSummary.readoutModels.signedDOG;
rectRM   = acrossOffsetSummary.readoutModels.rectifiedDOG;
hasSignedFit = isfield(signedRM, 'fit') && ~isempty(signedRM.fit) && ...
  isfield(signedRM.fit, 'fitSuccess') && signedRM.fit.fitSuccess;
if hasSignedFit
  hSigned = plot(signedRM.plotOffsetsDeg, signedRM.plotPredictedScale, 'm-', 'LineWidth', 1.2);
  DOGFitText(0.35, 0.98, 'Signed DOG', signedRM);
  legendHandles(end+1) = hSigned;
  legendLabels{end+1} = 'Fitted Signed DOG';
end
hasRectFit = isfield(rectRM, 'fit') && ~isempty(rectRM.fit) && ...
  isfield(rectRM.fit, 'fitSuccess') && rectRM.fit.fitSuccess;
if hasRectFit
  hRect = plot(rectRM.plotOffsetsDeg, rectRM.plotPredictedScale, 'b-', 'LineWidth', 1.2);
  DOGFitText(0.60, 0.98, 'Rectified DOG', rectRM);
  legendHandles(end+1) = hRect;
  legendLabels{end+1} = 'Fitted Rectified DOG';
end
legend(legendHandles, legendLabels, 'Location', 'northeast');
if hasSignedFit || hasRectFit
    title(sprintf('DOG Fits to Normalized Scales (%d bootstraps)', opts.NBoot));
else
    title(sprintf('Normalized Scales (No Fit Over %d bootstraps)', opts.NBoot));
end
scaleText(0.98, 0.88, offsets, obsScale, ci95, emp);
xlabel('Probe Offset (deg)');
ylabel('Normalized Scale');
xlim([0, 180]);
box off;
saveas(fig1, fullfile(opts.PlotDir, 'ScaleFits.pdf'));

% ---- Plots 2/3: fitted readout over MT preferred direction ----

if hasSignedFit
    tmpSummary = acrossOffsetSummary;
    tmpSummary.readoutModel = signedRM;
    plotReadoutDiagnostics(301, tmpSummary, opts);
end
if hasRectFit
  tmpSummary = acrossOffsetSummary;
  tmpSummary.readoutModel = rectRM;
  plotReadoutDiagnostics(302, tmpSummary, opts);
end
end

% ========================================================================
function DOGFitText(x, y, label, rm)
% One-model parameter/goodness-of-fit text block.

lines = {label};
for p = 1:min(numel(rm.params), numel(rm.paramNames))
    lines{end+1} = sprintf('  %s = %.4g', rm.paramNames{p}, rm.params(p)); %#ok<AGROW>
end
if isfield(rm, 'fit') && ~isempty(rm.fit) && isfield(rm.fit, 'goodnessOfFit') && isstruct(rm.fit.goodnessOfFit)
    g = rm.fit.goodnessOfFit;
    if isfield(g, 'weightedLoss') && isfinite(g.weightedLoss)
        lines{end+1} = sprintf('  loss = %.4g', g.weightedLoss); 
    end
    if isfield(g, 'reducedChiSq') && isfinite(g.reducedChiSq)
      lines{end+1} = sprintf('  red chi2 = %.4g', g.reducedChiSq); 
    end
    if isfield(g, 'aicc') && isfinite(g.aicc)
        lines{end+1} = sprintf('  AICc = %.4g', g.aicc); 
    end
end
txt = strjoin(lines, newline);
text(x, y, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 8, ...
  'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4, 'Interpreter', 'none');
end

% ========================================================================
function y = valueOrNaN(x)
if isempty(x)
  y = [NaN NaN];
else
  y = double(x);
end
end

% ========================================================================
function scaleText(x, y, offsets, obsScale, ci95, emp)
% One-model parameter/goodness-of-fit text block.

C = arrayfun(@(x) valueOrNaN(x.nTrialsByStep), emp, 'UniformOutput', false);
nTrialsByStep = vertcat(C{:});
stepTypes = [emp.stepType];
lines = {};
for index = 1:numel(offsets)
  if isnan(obsScale(index))
    continue;
  end
  nTrials = nTrialsByStep(index, stepTypes(index));
  lines{end+1} = sprintf('%3d°: scale %.2f, %.2f-%.2f 95%% CI (n = %5d)', ...
    offsets(index), obsScale(index), ci95(index, 1),  ci95(index, 2), nTrials); %#ok<AGROW>
end
txt = strjoin(lines, newline);
text(x, y, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', 'FontSize', 8, ...
  'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7], 'Margin', 4, 'Interpreter', 'none');
end

% ========================================================================
function val = getFieldOrDefault(S, fieldName, defaultVal)

if isfield(S, fieldName)
  val = S.(fieldName);
else
  val = defaultVal;
end
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

assert(isfield(sessionProbeHeader, 'probeDirDeg'), ...
    'updateAcrossOffsetSummaries:MissingProbeDir', ...
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

% ========================================================================
function files = pathsToDirStruct(filePaths)
% pathsToDirStruct  Convert full path cell/string array to dir-like struct array.
%
% updateAcrossOffsetSummaries historically iterates over dir() output with
% .folder and .name fields. selectAnalysisFiles returns selected full paths,
% so this adaptor preserves the existing downstream loop.

if isempty(filePaths)
    files = struct('folder', {}, 'name', {});
    return;
end

if isstring(filePaths)
    filePaths = cellstr(filePaths);
end

files = repmat(struct('folder', '', 'name', ''), numel(filePaths), 1);

for i = 1:numel(filePaths)
    [folderName, baseName, ext] = fileparts(char(filePaths{i}));
    files(i).folder = folderName;
    files(i).name = [baseName ext];
end
end