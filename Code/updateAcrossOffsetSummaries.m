function acrossOffsetSummary = updateAcrossOffsetSummaries(summaryDir, varargin)
% REQUIRED SESSION-LEVEL CONTRACT
%   This function operates on per-session kernel summary files produced by
%   compileKernelSessionSummary.
%
%   Each summary should contain:
%       sessionProbeHeader       authoritative metadata for this derived
%                                probe-offset session
%       probeOffsetDeg           numeric scalar probe direction / offset
%       sessionID/sessionName    session identifier
%       compStats                kernel comparison statistics
%       hitStats                 trial/hit counts
%       kernelFile               path to kernel file
%       noiseFile                path to noise-matrix file
%
%   The noiseFile should contain the resampling-ready fields accepted by
%   computeSessionKernels:
%       sessionProbeHeader
%       prefNoiseByPatch
%       probeNoiseByPatch
%       trialOutcomesAll
%       changeSidesAll
%       changeIndicesAll
%       lr
%       sideTypeNames
%
%   The original parent file header is not part of the downstream contract.
%   Any remaining support for header/blockStatus access is transitional
%   compatibility only and should not be used by newly generated derived
%   files.
%
% READOUT-FIT ELIGIBILITY
%   Kernel generation and averaging include all sessions not excluded by
%   excludeFile.m. This function applies additional readout-fit exclusions
%   internally, currently requiring paired-probe offsets:
%       0 < probeOffsetDeg < 180
%   and paired-probe noise amplitudes:
%       probeCohNoisePC == 7 or 10/sqrt(2), within tolerance.
%
%   Thus legacy 180 deg single-stream sessions may be processed and averaged
%   upstream, but are excluded here from the paired-probe readout fit.
% This function loads previously saved per-session summary files, groups them
% by probe offset, applies readout-fit eligibility criteria, performs
% across-offset bootstrap resampling, fits the readout model, and saves a
% single across-offset summary structure.
% -------------------------------------------------------------------------
    
if nargin < 1 || isempty(summaryDir)
  summaryDir = fullfile(folderPath(), 'Data');
end
opts = parseInputs(summaryDir, varargin{:});

if ~isempty(opts.RandomSeed)
    rng(opts.RandomSeed);
end

if opts.Verbose
    fprintf('updateAcrossOffsetSummaries: scanning %s\n', summaryDir);
end

sessionList = loadSessionSummaries(opts);
acrossOffsetSummary = initializeAcrossOffsetSummary(opts, sessionList);
if isempty(sessionList)
    warning('No usable session summaries found in %s.', summaryDir);
    saveAcrossOffsetSummary(opts, acrossOffsetSummary);
    return;
end

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
empirical = computeEmpiricalSummaries(offsetData, opts);
acrossOffsetSummary.empirical = empirical;

bootstrap = runHierarchicalBootstrap(offsetData, opts);
acrossOffsetSummary.bootstrap = bootstrap;

%% ---- Effective MT-to-choice weighting fit ----
offsetsDegAll   = [empirical.probeOffsetDeg];
pooledScaleAll  = [empirical.pooledScale];
% Bootstrap variance of the pooled scale at each measured offset. These are
% used as inverse-variance weights in fitReadoutDOGToScales.
bootstrapVarAll = bootstrap.fitBootstrap.offsetFitVar;
assert(numel(bootstrapVarAll) == numel(offsetsDegAll), 'bootstrap offset variances do not match empirical offsets');

% Fixed MT forward model used to map the DOG readout onto predicted scale.
mtModel = makeMTReadoutForwardModel('sigmaMTDeg', 37.5, 'phiDeg', -180:1:179);

% Store measured pooled scales by offset for downstream plotting/fits.
acrossOffsetSummary.measurements = struct();
acrossOffsetSummary.measurements.offsetsDeg   = offsetsDegAll(:)';
acrossOffsetSummary.measurements.pooledScale  = pooledScaleAll(:)';
acrossOffsetSummary.measurements.bootstrapVar = bootstrapVarAll(:)';

% Fit only non-anchor offsets. The 0-deg value is the normalization anchor.
isAnchor = abs(offsetsDegAll) < 1e-9;
fitOffsetsDeg = offsetsDegAll(~isAnchor);
fitScales     = pooledScaleAll(~isAnchor);
fitVars       = bootstrapVarAll(~isAnchor);

rm = struct();
rm.activeModelName = 'dog_readout';
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

if numel(fitOffsetsDeg) >= 1 && ...
        all(isfinite(fitScales)) && ...
        all(isfinite(fitVars)) && ...
        all(fitVars > 0)

    readoutFit = fitReadoutDOGToScales( ...
        fitOffsetsDeg, fitScales, fitVars, mtModel, ...
        'Bounds', opts.Bounds);

    rm.fit = readoutFit;
    rm.nFreeParams = readoutFit.nFreeParams;
    rm.phiDeg = mtModel.phiDeg;
    rm.readoutPhiRaw = readoutFit.readoutPhiRaw;
    rm.readoutPhi = readoutFit.readoutPhi;
    rm.paramNames = readoutFit.paramNames;
    rm.params = readoutFit.params;
    rm.paramStruct = readoutFit.paramStruct;
    rm.predictedAtMeasuredOffsets = ...
        predictNormalizedScaleFromReadout(readoutFit.params, offsetsDegAll, mtModel);
    rm.plotOffsetsDeg = 0:1:180;
    rm.plotPredictedScale = ...
        predictNormalizedScaleFromReadout(readoutFit.params, rm.plotOffsetsDeg, mtModel);

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
acrossOffsetSummary.readoutModel = rm;
acrossOffsetSummary.history = updateHistory(acrossOffsetSummary, opts);

saveAcrossOffsetSummary(opts, acrossOffsetSummary);
if opts.Verbose
    fprintf('updateAcrossOffsetSummaries: done. Saved summary to %s\n', opts.SaveFile);
end

if opts.MakePlots
    try
        makeAcrossOffsetPlots(acrossOffsetSummary, opts);
    catch ME
        warning('Plot generation failed: %s', ME.message);
    end
end

end

% ========================================================================
function opts = parseInputs(summaryDir, varargin)

if iscell(summaryDir)
    assert(~isempty(summaryDir), 'summaryDir cell array cannot be empty.');
elseif isstring(summaryDir)
    summaryDir = cellstr(summaryDir);
elseif ischar(summaryDir)
    % keep as char for now
else
    error('summaryDir must be a char, string, or cell array of paths.');
end
p = inputParser;
p.FunctionName = mfilename;


addRequired(p, 'summaryDir', @(x) ischar(x) || isstring(x) || iscell(x));
defaultSaveFile = fullfile(folderPath(), 'Data', 'AcrossOffsetSummaries', 'IDR_acrossOffsetSummary.mat');
addParameter(p, 'SaveFile', defaultSaveFile, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotDir',  fullfile(dataFolderPath(), '..', 'Plots', 'ReadoutFits'), @(x) ischar(x) || isstring(x));
addParameter(p, 'NBoot', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CILevels', [68 95], @(x) isnumeric(x) && isvector(x) && all(x > 0) && all(x < 100));
addParameter(p, 'Verbose', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MakePlots', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FilePattern', '*.mat', @(x) ischar(x) || isstring(x));
addParameter(p, 'OffsetField', 'probeOffsetDeg', @(x) ischar(x) || isstring(x));
addParameter(p, 'SessionNameField', 'sessionName', @(x) ischar(x) || isstring(x));
addParameter(p, 'ScaleMode', 'scaleFit', @(x) ischar(x) || isstring(x));
addParameter(p, 'ScaleSideType', 'change', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'ScaleStepType', 'inc', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'Bounds', struct(), @(x) isstruct(x));
addParameter(p, 'RandomSeed', [], @(x) isempty(x) || (isscalar(x) && isnumeric(x)));

parse(p, summaryDir, varargin{:});
opts = p.Results;
opts.summaryDir = normalizeSummaryDirs(summaryDir);
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

end

% ========================================================================
function sessionList = loadSessionSummaries(opts)
% Load per-session summaries and apply exclusion logic.

files = [];
for iDir = 1:numel(opts.summaryDir)
    thisDir = opts.summaryDir{iDir};
    files = [files; dir(fullfile(thisDir, opts.FilePattern))]; %#ok<AGROW>
end

sessionList = repmat(emptySessionStruct(), 0, 1);

if opts.Verbose
    fprintf('  Found %d candidate session summary files across %d folders.\n', ...
        numel(files), numel(opts.summaryDir));
end

for iFile = 1:numel(files)
    thisFile = fullfile(files(iFile).folder, files(iFile).name);
    try
        S = load(thisFile);
    catch ME
        warning('Could not load %s: %s', thisFile, ME.message);
        continue;
    end
    sessionStruct = extractSessionStruct(S, thisFile, opts);
    if isempty(sessionStruct)
        continue;
    end
    [tfExclude, reasonStr] = excludeFromReadoutFit(sessionStruct);
    sessionStruct.isExcluded = tfExclude;
    sessionStruct.excludeReason = reasonStr;
    sessionList(end+1,1) = sessionStruct; %#ok<AGROW>
end

if opts.Verbose
    nExcluded = sum([sessionList.isExcluded]);
    fprintf('  Loaded %d session summaries (%d excluded, %d included).\n', ...
        numel(sessionList), nExcluded, numel(sessionList) - nExcluded);
    if ~isempty(sessionList)
        fprintf('  Session diagnostics:\n');
        for iS = 1:numel(sessionList)
            incTxt = 'IN';
            if sessionList(iS).isExcluded
                incTxt = 'OUT';
            end
            fprintf('    %-3s  %-24s  offset=%6.1f  nTrials=%5g', ...
                incTxt, sessionList(iS).sessionName, sessionList(iS).probeOffsetDeg, sessionList(iS).nTrialsTotal);
            if sessionList(iS).isExcluded && ~isempty(sessionList(iS).excludeReason)
                fprintf('  reason=%s', sessionList(iS).excludeReason);
            end
            fprintf('\n');
        end
    end
end

end

% ========================================================================
function sessionStruct = extractSessionStruct(S, sourceFile, opts)
% Attempt to find a plausible session summary struct in a loaded MAT file.

sessionStruct = [];
fn = fieldnames(S);

% Preferred: variable already named summary (compileKernelSessionSummary output).
if isfield(S, 'summary') && isstruct(S.summary)
    sessionStruct = standardizeSessionStruct(S.summary, sourceFile, opts);
    return;
end

% Backward compatibility: variable named sessionSummary.
if isfield(S, 'sessionSummary') && isstruct(S.sessionSummary)
    sessionStruct = standardizeSessionStruct(S.sessionSummary, sourceFile, opts);
    return;
end

% Otherwise find the first struct containing plausible fields.
for i = 1:numel(fn)
    val = S.(fn{i});
    if isstruct(val)
        if isfield(val, 'scale') || isfield(val, opts.OffsetField) || isfield(val, 'header')
            sessionStruct = standardizeSessionStruct(val, sourceFile, opts);
            return;
        end
    end
end

warning('No recognizable session summary struct found in %s.', sourceFile);

end

% ========================================================================
function sessionStruct = standardizeSessionStruct(inStruct, sourceFile, opts)
% Map session summary into a stable internal struct.

sessionStruct = emptySessionStruct();

sessionStruct.probeOffsetDeg = extractProbeOffsetDeg(inStruct, opts.OffsetField);
sessionStruct.sessionName    = extractSessionName(inStruct, sourceFile, opts.SessionNameField);
sessionStruct.sessionDate    = extractSessionDate(inStruct, sourceFile);
sessionStruct.header = getFieldOrDefault(inStruct, 'sessionProbeHeader', struct());

if isempty(fieldnames(sessionStruct.header))
    sessionStruct.header = getFieldOrDefault(inStruct, 'header', struct());
end
if isempty(fieldnames(sessionStruct.header)) && ...
        isfield(inStruct, 'kernelFile') && ...
        (ischar(inStruct.kernelFile) || isstring(inStruct.kernelFile))
    try
        K = load(char(inStruct.kernelFile), 'sessionProbeHeader');
        if isfield(K, 'sessionProbeHeader') && isstruct(K.sessionProbeHeader)
            sessionStruct.header = K.sessionProbeHeader;
        end
    catch
    end
end
sessionStruct.compStats      = getFieldOrDefault(inStruct, 'compStats', struct());
sessionStruct.hitStats       = getFieldOrDefault(inStruct, 'hitStats', struct());
sessionStruct.bootstrapSource = extractBootstrapSource(inStruct);
sessionStruct.rawSummary     = inStruct;
sessionStruct.scalePointEstimate = extractSessionScaleEstimate(inStruct, opts);
try
    sessionStruct.scalePointEstimate = recomputeSessionScalePointEstimate(sessionStruct, opts);
catch
end
sessionStruct.nTrialsTotal = extractSessionTrialCount(inStruct);
if isfield(inStruct, 'flags') && isstruct(inStruct.flags)
    if isfield(inStruct.flags, 'excluded') && islogical(inStruct.flags.excluded)
        sessionStruct.isExcluded = inStruct.flags.excluded;
        if sessionStruct.isExcluded
            sessionStruct.excludeReason = 'summary.flags.excluded';
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
    'header', struct(), ...
    'compStats', struct(), ...
    'hitStats', struct(), ...
    'bootstrapSource', struct(), ...
    'rawSummary', struct(), ...
    'scalePointEstimate', NaN, ...
    'nTrialsTotal', NaN, ...
    'isExcluded', false, ...
    'excludeReason', '', ...
    'sourceFile', '' );
end

% ========================================================================
function [exclude, reason] = excludeFromReadoutFit(sessionStruct)
% excludeFromReadoutFit  Exclude sessions not eligible for paired-probe readout fitting.
%
% These exclusions are specific to the across-offset readout model. They do
% not imply that the session should be excluded from kernel generation,
% kernel averaging, or other preprocessing summaries.

exclude = false;
reason = '';

H = sessionStruct.header;

% Prefer the standardized offset field, because standardizeSessionStruct()
% has already extracted it from the summary metadata.
probeDirDeg = sessionStruct.probeOffsetDeg;

if ~(isfinite(probeDirDeg) && probeDirDeg > 0 && probeDirDeg < 180)
    exclude = true;
    reason = sprintf('probeDirDeg %.6g is not a paired-probe offset', probeDirDeg);
    return;
end

[probeCohNoisePC, ok] = localGetHeaderValue(H, 'probeCohNoisePC');
if ~ok || ~isfinite(probeCohNoisePC)
    exclude = true;
    reason = 'missing/non-finite probeCohNoisePC';
    return;
end

targetProbeNoisePCs = [7, 10 / sqrt(2)];
tol = 1e-4;

if all(abs(probeCohNoisePC - targetProbeNoisePCs) > tol)
    exclude = true;
    reason = sprintf('probeCohNoisePC %.6g is not paired-probe amplitude', probeCohNoisePC);
    return;
end
end

% ========================================================================
function [v, ok] = localGetHeaderValue(S, fieldName)

ok = false;
v = NaN;

if isfield(S, fieldName)
    x = S.(fieldName);
    if isstruct(x) && isfield(x, 'data')
        v = double(x.data(1));
    else
        v = double(x(1));
    end
    ok = true;
    return;
end

if isfield(S, 'blockStatus') && isfield(S.blockStatus, 'data') && ...
        isfield(S.blockStatus.data, fieldName)

    x = S.blockStatus.data(1).(fieldName);
    v = double(localExtractValue(x));
    ok = true;
    return;
end
end

% ========================================================================
function acrossOffsetSummary = initializeAcrossOffsetSummary(opts, sessionList)

acrossOffsetSummary = struct();
acrossOffsetSummary.meta = struct( ...
    'analysisDate', datetime('now'), ...
    'summaryDirs', {opts.summaryDir}, ...
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
acrossOffsetSummary.modelFitsNote = ...
    ['Legacy offset-space modelFits removed. Primary interpretation uses ' ...
     'acrossOffsetSummary.readoutModel.'];
end

% ========================================================================
function [offsetData, offsetKeys, usedFiles] = buildOffsetData(sessionList, ~)
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
    offsetData(k).nTrialsBySession = [included.nTrialsTotal];
    offsetData(k).nTrialsTotal = sum(offsetData(k).nTrialsBySession, 'omitnan');
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
    'nTrialsTotal', 0, ...
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
    'pooledScale', NaN, ...
    'meanSessionScale', NaN, ...
    'medianSessionScale', NaN, ...
    'semSessionScale', NaN, ...
    'boot68', [NaN NaN], ...
    'boot95', [NaN NaN], ...
    'nSessions', 0, ...
    'nTrialsTotal', 0, ...
    'sessionWeighting', 'inverse_variance_kernel_pool', ...
    'distributionSkew', NaN), numel(offsetData), 1);

for k = 1:numel(offsetData)
    x = offsetData(k).scaleBySession;
    empirical(k).probeOffsetDeg   = offsetData(k).probeOffsetDeg;
    empirical(k).pooledScale      = computeOffsetPooledScale(offsetData(k), opts, false);
    empirical(k).meanSessionScale = mean(x, 'omitnan');
    empirical(k).medianSessionScale = median(x, 'omitnan');
    empirical(k).semSessionScale  = localSEM(x);
    empirical(k).boot68           = percentileCI(x, 68);
    empirical(k).boot95           = percentileCI(x, 95);
    empirical(k).nSessions        = offsetData(k).nSessions;
    empirical(k).nTrialsTotal     = offsetData(k).nTrialsTotal;
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
end

end

% ========================================================================
function bootstrap = runHierarchicalBootstrap(offsetData, opts)
% Hierarchical bootstrap over sessions and trials.
%
% Purpose:
%   Estimate uncertainty in the pooled psychophysical scale at each tested
%   probe offset. These offset variances are then used as weights for the
%   primary MT-readout DOG fit.
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
             'the primary readout DOG fit.'] );
end

% ========================================================================
function scaleVal = recomputeSessionScalePointEstimate(sessionStruct, opts)
% Recompute the requested side/step point estimate from source data.

bs = sessionStruct.bootstrapSource;
if isempty(fieldnames(bs))
    scaleVal = sessionStruct.scalePointEstimate;
    return;
end

if isfield(bs, 'sessionData') && isstruct(bs.sessionData)
    sessionData = bs.sessionData;
elseif hasComputeSessionKernelFields(bs)
    sessionData = bs;
elseif isfield(bs, 'noiseFile') && exist(bs.noiseFile, 'file')
    sessionData = load(bs.noiseFile);
else
    scaleVal = sessionStruct.scalePointEstimate;
    return;
end

[~, ~, ~, ~, compStats] = computeSessionKernels(sessionData, []);

if isfield(compStats, 'normScale')
    scaleVal = selectCompStatsEntry(compStats.normScale, opts.ScaleSideType, opts.ScaleStepType);
else
    warning('updateAcrossOffsetSummaries:MissingNormScale', ...
        'compStats.normScale missing for session %s; falling back to raw scale.', ...
        sessionStruct.sessionName);
    scaleVal = selectCompStatsEntry(compStats.scale, opts.ScaleSideType, opts.ScaleStepType);
end
end

% ========================================================================
function [kernels, kVars, hitStats, compStats] = recomputeSessionKernelStruct(sessionStruct, ~, doTrialBootstrap)
% Recompute full session kernel outputs from source data.

bs = sessionStruct.bootstrapSource;

if isfield(bs, 'sessionData') && isstruct(bs.sessionData)
    sessionData = bs.sessionData;
elseif hasComputeSessionKernelFields(bs)
    sessionData = bs;
elseif isfield(bs, 'noiseFile') && exist(bs.noiseFile, 'file')
    sessionData = load(bs.noiseFile);
else
    error('Missing bootstrap source for session %s.', sessionStruct.sessionName);
end

trialIdx = [];
if doTrialBootstrap
    nTrials = size(sessionData.prefNoiseByPatch, 3);
    trialIdx = randi(nTrials, [1 nTrials]);
end
[kernels, kVars, ~, hitStats, compStats] = computeSessionKernels(sessionData, trialIdx);
end

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
    deltaM = (1 / sqrt(2)) .* (deltaMPlus + deltaMMinus);
end
end

% ========================================================================
function summaryDirs = normalizeSummaryDirs(summaryDir)

if ischar(summaryDir) || isstring(summaryDir)
    summaryDir = cellstr(summaryDir);
end

summaryDirs = {};
for i = 1:numel(summaryDir)
    d = char(summaryDir{i});

    % If the user passed Data/, discover Data/probe*/KernelSummaries.
    probeDirs = dir(fullfile(d, 'probe*'));
    probeDirs = probeDirs([probeDirs.isdir]);

    if ~isempty(probeDirs)
        for p = 1:numel(probeDirs)
            candidate = fullfile(probeDirs(p).folder, probeDirs(p).name, 'KernelSummaries');
            if exist(candidate, 'dir')
                summaryDirs{end+1} = candidate; %#ok<AGROW>
            end
        end

    % If the user passed a probe folder, append KernelSummaries.
    elseif exist(fullfile(d, 'KernelSummaries'), 'dir')
        summaryDirs{end+1} = fullfile(d, 'KernelSummaries'); %#ok<AGROW>

    % If the user passed KernelSummaries directly, use it.
    elseif exist(d, 'dir')
        summaryDirs{end+1} = d; %#ok<AGROW>

    else
        warning('updateAcrossOffsetSummaries:MissingSummaryDir', ...
            'Summary directory not found: %s', d);
    end
end

summaryDirs = unique(summaryDirs, 'stable');

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
    [kernels, kVars, ~, ~] = recomputeSessionKernelStruct(src, opts, doBootstrap);

    [kernelsNorm, kVarsNorm] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, src.header);

    sessionKernelsNorm{j} = kernelsNorm;
    sessionKVarsNorm{j}   = kVarsNorm;

    thisFrameRateHz = src.header.frameRateHz.data(1);
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
function tf = hasComputeSessionKernelFields(S)

req = {'prefNoiseByPatch','probeNoiseByPatch','trialOutcomesAll', ...
       'changeSidesAll','changeIndicesAll','sessionProbeHeader'};
tf = all(isfield(S, req));
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
% function diag = diagnoseReadoutDOGReachabilityN(offsetsDeg, targetScale, fitVars, varargin)
% 
% p = inputParser;
% addParameter(p, 'SigmaCenterDeg', [0.001 0.01 0.1 0.3 1 2 5 10 20 40 80 120], @(x) isnumeric(x) && isvector(x));
% addParameter(p, 'SigmaSurroundDeg', 1:2:300, @(x) isnumeric(x) && isvector(x));
% addParameter(p, 'SurroundGain', 0:0.01:2.0, @(x) isnumeric(x) && isvector(x));
% addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
% addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x));
% parse(p, varargin{:});
% opts = p.Results;
% 
% offsetsDeg = offsetsDeg(:)';
% targetScale = targetScale(:)';
% 
% mtModel = makeMTReadoutForwardModel( ...
%     'sigmaMTDeg', opts.SigmaMTDeg, ...
%     'phiDeg', opts.PhiDeg);
% 
% sC = opts.SigmaCenterDeg(:)';
% sS = opts.SigmaSurroundDeg(:)';
% aS = opts.SurroundGain(:)';
% 
% bestLoss = Inf;
% bestParams = [];
% bestPred = [];
% 
% for i = 1:numel(sC)
%     for j = 1:numel(sS)
%         if sS(j) < 1.25 * sC(i)
%             continue;
%         end
% 
%         for k = 1:numel(aS)
%             params = [sC(i), sS(j), aS(k)];
%             pred = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel);
% 
%             if ~all(isfinite(pred))
%                 continue;
%             end
% 
%             loss = sum(((pred - targetScale).^2) ./ fitVars);
%             if loss < bestLoss
%                 bestLoss = loss;
%                 bestParams = params;
%                 bestPred = pred;
%             end
%         end
%     end
% end
% 
% diag = struct();
% diag.offsetsDeg = offsetsDeg;
% diag.targetScale = targetScale;
% diag.bestParams = bestParams;
% diag.bestPrediction = bestPred;
% diag.bestLoss = bestLoss;
% diag.bestRMSE = sqrt(bestLoss / numel(targetScale));
% 
% fprintf('\nDOG reachability diagnostic:\n');
% disp(table(offsetsDeg(:), targetScale(:), bestPred(:), ...
%     'VariableNames', {'offsetDeg','target','bestDOG'}));
% fprintf('best params: sigmaC %.4g, sigmaS %.4g, As %.4g\n', bestParams);
% fprintf('best RMSE: %.4g\n', diag.bestRMSE);
% end

% ========================================================================
% function predScale = predictReadoutDOGScale(params, offsetsDeg, mtModel, baselineMode, fixedBaseline)
% 
% if nargin < 4 || isempty(baselineMode)
%     baselineMode = 'fixed';
% end
% 
% if nargin < 5 || isempty(fixedBaseline) || ~isfinite(fixedBaseline)
%     fixedBaseline = 0;
% end
% 
% baselineMode = lower(char(string(baselineMode)));
% 
% dogParams = params(1:3);
% 
% if strcmp(baselineMode, 'fit')
%     baseline = params(4);
% else
%     baseline = fixedBaseline;
% end
% 
% predScale = predictNormalizedScaleFromReadout(dogParams, offsetsDeg, mtModel) + baseline;
% end

% ========================================================================
function [aPhi, paramStruct] = evaluateReadoutDOG(phiDeg, params)
% Evaluate three-parameter DOG readout over MT preferred direction.
%
% Readout:
%   a(phi) = exp(-(phi^2)/(2*sigmaC^2)) ...
%          - As * exp(-(phi^2)/(2*sigmaS^2))
%
% params:
%   [sigmaC, sigmaS, As]

phiDeg = phiDeg(:)';

if numel(params) ~= 3
    error('evaluateReadoutDOG requires params = [sigmaC, sigmaS, As].');
end

sigmaC = params(1);
sigmaS = params(2);
As     = params(3);

aPhi = exp(-(phiDeg.^2) ./ (2 * sigmaC.^2)) ...
     - As .* exp(-(phiDeg.^2) ./ (2 * sigmaS.^2));

paramStruct = struct( ...
    'sigmaCenterDeg', sigmaC, ...
    'sigmaSurroundDeg', sigmaS, ...
    'surroundGain', As );
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
function predScale = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel)
% The probe template follows the stimulus construction:
%   - delta = 0 deg: single-channel template
%   - 0 < delta < 180 deg: two symmetric probe-noise streams at +delta and
%     -delta, each with amplitude 1/sqrt(2) relative to the single-stream
%     probe amplitude
%   - delta = 180 deg: single-channel template
%
% For 0 < delta < 180:
%   Delta m_eff(phi;delta) = (1/sqrt(2)) * [Delta m(phi;+delta) + Delta m(phi;-delta)]
%
% S_pred(delta) = sum_phi a(phi) Delta m_eff(phi;delta) / ...
%                 sum_phi a(phi) Delta m_eff(phi;0)

phiDeg = mtModel.phiDeg;
[aPhi, ~] = evaluateReadoutDOG(phiDeg, params);


refVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, 0, mtModel));

offsetsDeg = offsetsDeg(:)';
predScale = nan(size(offsetsDeg));
for i = 1:numel(offsetsDeg)
    probeVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, offsetsDeg(i), mtModel));
    if isfinite(refVal) && abs(refVal) > 0
        predScale(i) = probeVal / refVal;
    else
        predScale(i) = NaN;
    end
end
end

% ========================================================================
% function predScale = predictNormalizedScaleFromExplicitReadout(offsetsDeg, mtModel, aPhi)
% % Predict normalized scale from an explicit readout vector a(phi).
% %
% % The probe template follows the stimulus construction:
% %   - delta = 0 deg: single-channel template
% %   - 0 < delta < 180 deg: symmetrized template for equal probe-noise
% %     contributions at +delta and -delta
% %   - delta = 180 deg: single-channel template
% 
% phiDeg = mtModel.phiDeg(:)';
% aPhi   = aPhi(:)';
% 
% if numel(aPhi) ~= numel( ...
%     phiDeg)
%     error('predictNormalizedScaleFromExplicitReadout: aPhi must match mtModel.phiDeg in length.');
% end
% 
% refVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, 0, mtModel));
% tol = 1e-10;
% 
% offsetsDeg = offsetsDeg(:)';
% predScale = nan(size(offsetsDeg));
% 
% for i = 1:numel(offsetsDeg)
%     probeVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, offsetsDeg(i), mtModel));
%     if isfinite(refVal) && abs(refVal) > tol
%         predScale(i) = probeVal / refVal;
%     else
%         predScale(i) = NaN;
%     end
% end
% end

% ========================================================================
function plotReadoutDiagnostics(acrossOffsetSummary, opts)
% Plot fitted readout, MT templates, and their products to visualize how
% overlap determines predicted normalized scale.

  rm = acrossOffsetSummary.readoutModel;
  if ~isfield(rm, 'fit') || isempty(rm.fit) || ~rm.fit.fitSuccess
      return;
  end

  phiDeg = rm.phiDeg(:)';
  aPhi   = rm.readoutPhi(:)';   % normalized display readout, a(0)=1
  mtp = rm.mtForwardModelParams;
  mtModel = makeMTReadoutForwardModel('sigmaMTDeg', mtp.sigmaMTDeg, 'phiDeg', mtp.phiDeg);
  offsetsDeg = [0, [acrossOffsetSummary.offsetData.probeOffsetDeg]];
  nOffsets = numel(offsetsDeg);
  
  deltaM   = cell(1, nOffsets);
  prodTerm = cell(1, nOffsets);
  overlap  = nan(1, nOffsets);
  posPart  = nan(1, nOffsets);
  negPart  = nan(1, nOffsets);

  for i = 1:nOffsets
      deltaM{i} = computeMTSymmetrizedDeltaM(phiDeg, offsetsDeg(i), mtModel);
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
  fig = figure(301); clf;
  tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
  lineCol = lines(nOffsets);
  
  % ---- top panel: readout function ----
  nexttile; hold on;
  hFitReadout  = plot(phiDeg, aPhi, 'k-', 'LineWidth', 2);
  hFlatReadout = plot(phiDeg, ones(size(phiDeg)), 'k--', 'LineWidth', 1.5);
  plot(phiDeg, zeros(size(phiDeg)), 'k:', 'LineWidth', 1);
  xlabel('\phi (deg)');
  ylabel('a(\phi)');
  title(sprintf('Fitted DOG readout over MT preferred direction (%d bootstraps)', opts.NBoot));
  legend([hFitReadout, hFlatReadout], [{'Fitted readout a(\phi)'}, {'Flat readout = 1'}], 'Location', 'northeast');
  paramText = cell(rm.nFreeParams, 1);
  for p = 1:rm.nFreeParams
    paramText{p} = sprintf('%s: %.2f', rm.paramNames{p}, rm.params(p));
  end
  text(-120, 0.95, paramText, 'horizontalAlignment', 'right', 'VerticalAlignment', 'top');
  ylimits = ylim();
  ylim([min(0.2, ylimits(1)), max(1.1, ylimits(2))]);
  box off;
  
  % ---- middle panel: MT populations responses ----
  nexttile; hold on;
  title('MT Population Responses to Probes (Flat Readout)');
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
      saveas(fig, fullfile(opts.PlotDir, 'ReadoutFunctions.pdf'));
  end

  % Print a compact numeric summary to the command window
  % fprintf('\nReadout overlap diagnostic:\n');
  % for i = 1:nOffsets
  %     fprintf(['  delta = %6.1f deg:  <a,Delta m> = %9.4f   ' ...
  %              'positive = %9.4f   negative = %.2g'], ...
  %         offsetsDeg(i), overlap(i), posPart(i), negPart(i));
  %     if ~isnan(predScaleFit(i))
  %         fprintf('   S_fit(delta) = %.2g', predScaleFit(i));
  %     end
  %     if ~isnan(predScaleFlat(i))
  %         fprintf('   S_flat(delta) = %.2g', predScaleFlat(i));
  %     else
  %       fprintf('   S_flat(delta) = undefined (flat readout gives zero overlap)');
  %     end
  %     fprintf('\n');
  % end
end

% ========================================================================
function ang = wrapTo180Local(ang)
% Map angles to [-180, 180).

ang = mod(ang + 180, 360) - 180;
end

% ========================================================================
function fitResult = fitReadoutDOGToScales(offsetsDeg, obsScale, obsVar, mtModel, varargin)
% Fit three-parameter DOG readout to observed non-anchor scale values.
%
% Weighted objective:
%   sum_i (obsScale_i - predScale_i)^2 / obsVar_i

p = inputParser;
addParameter(p, 'Bounds', struct(), @(x) isstruct(x));
parse(p, varargin{:});
opts = p.Results;

[p0, lb, ub] = initialGuessReadoutDOG(offsetsDeg, opts.Bounds);

obj = @(params) readoutDOGObjective(params, offsetsDeg, obsScale, obsVar, mtModel);

p0 = p0(:).';
lb = lb(:).';
ub = ub(:).';
p0 = min(max(p0, lb), ub);

params = nan(size(p0));
loss = NaN;

try
    fminconOpts = optimoptions( ...
        'fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'interior-point');

    [params, loss, exitflag] = fmincon( ...
        obj, ...
        p0, ...
        [], [], [], [], ...
        lb, ub, ...
        [], ...
        fminconOpts);

    fitSuccess = exitflag > 0 && all(isfinite(params));
catch ME
    warning('fmincon failed in fitReadoutDOGToScales: %s', ME.message);
    fitSuccess = false;
end

predMeasured = nan(size(offsetsDeg));
paramStruct = struct();

if fitSuccess
    predMeasured = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel).';
    [~, paramStruct] = evaluateReadoutDOG(mtModel.phiDeg, params);
end

fitResult = struct();
fitResult.modelName = 'dog_readout';
fitResult.fitSuccess = fitSuccess;
fitResult.loss = loss;
fitResult.params = params;
fitResult.paramStruct = paramStruct;
fitResult.paramNames = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};
fitResult.nFreeParams = 3;
fitResult.offsetsDeg = offsetsDeg(:)';
fitResult.observedScale = obsScale(:)';
fitResult.observedVar = obsVar(:)';
fitResult.predictedScale = predMeasured(:)';
fitResult.phiDeg = mtModel.phiDeg;

if fitSuccess
    fitResult.readoutPhiRaw = evaluateReadoutDOG(mtModel.phiDeg, params);
    fitResult.readoutPhi = normalizeReadoutAtPreferred(mtModel.phiDeg, fitResult.readoutPhiRaw);
else
    fitResult.readoutPhiRaw = nan(size(mtModel.phiDeg));
    fitResult.readoutPhi = nan(size(mtModel.phiDeg));
end
end

% ========================================================================
function [p0, lb, ub] = initialGuessReadoutDOG(offsetsDeg, bounds)

sigmaC0 = max(10, min(50, median(offsetsDeg(offsetsDeg > 0), 'omitnan')));
if isempty(sigmaC0) || ~isfinite(sigmaC0)
    sigmaC0 = 25;
end

sigmaS0 = max(sigmaC0 + 20, 90);
As0 = 0.5;

p0 = [sigmaC0, sigmaS0, As0];
lb = [1e-3, 1e-3, 0];
ub = [300, 300, 10];

if isfield(bounds, 'readoutDOG')
    B = bounds.readoutDOG;
elseif isfield(bounds, 'dog_readout')
    B = bounds.dog_readout;
elseif isfield(bounds, 'dog')
    B = bounds.dog;
else
    B = struct();
end

if isfield(B, 'lb'), lb = B.lb; end
if isfield(B, 'ub'), ub = B.ub; end

p0 = p0(:);
lb = lb(:);
ub = ub(:);
end

% ========================================================================
function sse = readoutDOGObjective(params, offsetsDeg, obsScale, obsVar, mtModel)

predScale = predictNormalizedScaleFromReadout(params, offsetsDeg, mtModel);
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
% Primary plots for the DOG readout model:
%   1) observed scale values by probe offset, with predicted values overlaid
%      when a readout fit is available
%   2) fitted readout a(phi) over MT preferred direction, only when fit exists

emp = acrossOffsetSummary.empirical;
offsets = [emp.probeOffsetDeg];
obsScale = [emp.pooledScale];
ci68 = vertcat(emp.boot68);

hasReadoutFit = isfield(acrossOffsetSummary, 'readoutModel') && ...
    isfield(acrossOffsetSummary.readoutModel, 'fit') && ~isempty(acrossOffsetSummary.readoutModel.fit) && ...
    isfield(acrossOffsetSummary.readoutModel.fit, 'fitSuccess') && acrossOffsetSummary.readoutModel.fit.fitSuccess;

% ---- Plot 1: observed scale by offset; overlay prediction if available ----
fig1 = figure(300); clf; hold on;
hObs = errorbar([0, offsets], [1, obsScale], [0, obsScale - ci68(:,1)'], [0, ci68(:,2)' - obsScale], ...
    'bo-', 'LineWidth', 1.2, 'MarkerFaceColor', 'b');
% plot(0, 1, 'ko', 'MarkerFaceColor', 'k');
plot([0, 180], [0,0], 'k:');
if hasReadoutFit
    rm = acrossOffsetSummary.readoutModel;
    predScale = rm.predictedAtMeasuredOffsets;
    hPred = plot(offsets, predScale, 'k-', 'LineWidth', 2);
    legend([hObs, hPred], {'Observed pooled scale (68% CI)', 'Predicted from fitted readout'}, 'Location', 'northeast');
    title(sprintf('Fit to Normalized Scales (%d bootstraps)', opts.NBoot));
else
    legend(hObs, {'Observed Pooled Normalized Scale (68% CI)'}, 'Location', 'northeast');
    title(sprintf('Normalized Scales (No Fit Over %d bootstraps)', opts.NBoot));
end
xlabel('Probe Offset (deg)');
ylabel('Normalized Scale');
xlim([0, 180])
box off;
saveas(fig1, fullfile(opts.PlotDir, 'ScaleFits.pdf'));

% ---- Plot 2: fitted readout over MT preferred direction ----
if hasReadoutFit
  plotReadoutDiagnostics(acrossOffsetSummary, opts);
end

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
function name = extractSessionName(inStruct, sourceFile, sessionNameField)

name = '';

if isfield(inStruct, sessionNameField)
    try
        name = char(string(inStruct.(sessionNameField)));
    catch
    end
end

if isempty(name) && isfield(inStruct, 'sessionID')
    name = char(string(inStruct.sessionID));
end

if isempty(name)
    [~, name, ~] = fileparts(sourceFile);
end

end

% ========================================================================
function scaleVal = extractSessionScaleEstimate(inStruct, opts)
% Project-specific adaptor. Tries common fields and then falls back.

scaleVal = NaN;

if isfield(inStruct, 'scalePointEstimate') && isfinite(inStruct.scalePointEstimate)
    scaleVal = inStruct.scalePointEstimate;
    return;
end
if isfield(inStruct, 'scale') && isstruct(inStruct.scale)
    if isfield(inStruct.scale, 'normEstimate') && isfinite(inStruct.scale.normEstimate)
        scaleVal = inStruct.scale.normEstimate;
        return;
    end
    if isfield(inStruct.scale, 'estimate') && isfinite(inStruct.scale.estimate)
        scaleVal = inStruct.scale.estimate;
        return;
    end
end
if isfield(inStruct, 'summary') && isstruct(inStruct.summary)
    if isfield(inStruct.summary, 'scalePointEstimate') && isfinite(inStruct.summary.scalePointEstimate)
        scaleVal = inStruct.summary.scalePointEstimate;
        return;
    end
end
if isfield(inStruct, 'compStats') && isstruct(inStruct.compStats)
    cs = inStruct.compStats;
    if isfield(cs, 'normScale') && isnumeric(cs.normScale)
        try
            scaleVal = selectCompStatsEntry(cs.normScale, opts.ScaleSideType, opts.ScaleStepType);
            return;
        catch
        end
    end
    if isfield(cs, 'scale') && isnumeric(cs.scale)
        try
            scaleVal = selectCompStatsEntry(cs.scale, opts.ScaleSideType, opts.ScaleStepType);
            return;
        catch
        end
    end
end
if isfield(inStruct, 'avgCompStats') && isstruct(inStruct.avgCompStats)
    cs = inStruct.avgCompStats;
    if isfield(cs, 'normScale') && isnumeric(cs.normScale)
        try
            scaleVal = selectCompStatsEntry(cs.normScale, opts.ScaleSideType, opts.ScaleStepType);
            return;
        catch
        end
    end
    if isfield(cs, 'scale') && isnumeric(cs.scale)
        try
            scaleVal = selectCompStatsEntry(cs.scale, opts.ScaleSideType, opts.ScaleStepType);
            return;
        catch
        end
    end
end
end

% ========================================================================
function nTrials = extractSessionTrialCount(inStruct)

nTrials = NaN;

if isfield(inStruct, 'metrics') && isstruct(inStruct.metrics)
    if isfield(inStruct.metrics, 'nTrials') && isnumeric(inStruct.metrics.nTrials)
        nTrials = inStruct.metrics.nTrials;
        return;
    end
end

if isfield(inStruct, 'nTrialsTotal') && isnumeric(inStruct.nTrialsTotal)
    nTrials = inStruct.nTrialsTotal;
    return;
end

if isfield(inStruct, 'avgHitStats') && isstruct(inStruct.avgHitStats)
    hs = inStruct.avgHitStats;
    if isfield(hs, 'nTrials') && isnumeric(hs.nTrials)
        nTrials = sum(hs.nTrials(:), 'omitnan');
        return;
    end
end

if isfield(inStruct, 'hitStats') && isstruct(inStruct.hitStats)
    hs = inStruct.hitStats;
    if isfield(hs, 'nTrials') && isnumeric(hs.nTrials)
        nTrials = sum(hs.nTrials(:), 'omitnan');
        return;
    end
end

end

% ========================================================================
function dt = extractSessionDate(inStruct, sourceFile)

dt = NaT;

if isfield(inStruct, 'date')
    try
        dt = datetime(inStruct.date);
        return;
    catch
    end
end

if isfield(inStruct, 'sessionDate')
    try
        dt = datetime(inStruct.sessionDate);
        return;
    catch
    end
end

% try parsing YYYYMMDD from filename
expr = '(20\d{6})';
tok = regexp(sourceFile, expr, 'tokens', 'once');
if ~isempty(tok)
    try
        dt = datetime(tok{1}, 'InputFormat', 'yyyyMMdd');
    catch
    end
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
function probeOffsetDeg = extractProbeOffsetDeg(inStruct, offsetField)

probeOffsetDeg = NaN;

% New preferred contract: summary.sessionProbeHeader.probeDirDeg
if isfield(inStruct, 'sessionProbeHeader') && isstruct(inStruct.sessionProbeHeader) && ...
        isfield(inStruct.sessionProbeHeader, 'probeDirDeg')
    probeOffsetDeg = double(localExtractValue(inStruct.sessionProbeHeader.probeDirDeg));
    assert(isscalar(probeOffsetDeg) && isfinite(probeOffsetDeg), ...
        'sessionProbeHeader.probeDirDeg must be a finite scalar.');
    return;
end

% New explicit summary-level field, if present.
if isfield(inStruct, offsetField)
    probeOffsetDeg = double(localExtractValue(inStruct.(offsetField)));
    assert(isscalar(probeOffsetDeg) && isfinite(probeOffsetDeg), ...
        '%s must be a finite scalar.', offsetField);
    return;
end

% Transitional compatibility: summary.header may now actually be a
% sessionProbeHeader-like struct.
if isfield(inStruct, 'header') && isstruct(inStruct.header) && ...
        isfield(inStruct.header, 'probeDirDeg')
    probeOffsetDeg = double(localExtractValue(inStruct.header.probeDirDeg));
    assert(isscalar(probeOffsetDeg) && isfinite(probeOffsetDeg), ...
        'header.probeDirDeg must be a finite scalar.');
    return;
end

error('extractProbeOffsetDeg:MissingProbeOffset', ...
    'No authoritative probe offset found in summary sessionProbeHeader, header, or %s.', offsetField);
end

% ========================================================================
function val = localExtractValue(x)
% localExtractValue  Pull scalar/string value from either scalar or .data field.

if isstruct(x) && isfield(x, 'data')
    val = x.data;
else
    val = x;
end

if isnumeric(val) && ~isscalar(val)
    val = val(1);
end
end

% ========================================================================
function bs = extractBootstrapSource(inStruct)

bs = struct();

candidateFields = {'bootstrapSource','noiseMatrices','sessionKernelSource','resampleSource'};
for i = 1:numel(candidateFields)
    if isfield(inStruct, candidateFields{i}) && isstruct(inStruct.(candidateFields{i}))
        bs = inStruct.(candidateFields{i});
        return;
    end
end

% compileKernelSessionSummary stores the path to the noise matrix file.
if isfield(inStruct, 'noiseFile') && (ischar(inStruct.noiseFile) || isstring(inStruct.noiseFile))
    bs.noiseFile = char(inStruct.noiseFile);
end
if isfield(inStruct, 'kernelFile') && (ischar(inStruct.kernelFile) || isstring(inStruct.kernelFile))
    bs.kernelFile = char(inStruct.kernelFile);
end
if isfield(inStruct, 'track') && isstruct(inStruct.track)
    bs.track = inStruct.track;
end
if isfield(inStruct, 'header') && isstruct(inStruct.header)
    bs.header = inStruct.header;
end

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
function [kernelsNorm, kVarsNorm, normInfo] = normalizeProbeKernelsToPrefAmplitude(kernels, kVars, sessionProbeHeader)

normInfo = normalizationInfoFromSessionProbeHeader(sessionProbeHeader);
probeNormFactor = normInfo.probeNormFactor;

kernelsNorm = kernels;
kVarsNorm   = kVars;

kernelsNorm(:, :, 2, :) = kernelsNorm(:, :, 2, :) * probeNormFactor;
kVarsNorm(:, :, 2)      = kVarsNorm(:, :, 2) * probeNormFactor^2;
end

% ========================================================================
function normInfo = normalizationInfoFromSessionProbeHeader(sessionProbeHeader)

[prefCohNoisePC, okPref] = localGetHeaderValue(sessionProbeHeader, 'prefCohNoisePC');
[probeCohNoisePC, okProbe] = localGetHeaderValue(sessionProbeHeader, 'probeCohNoisePC');

if ~okPref || ~okProbe || ...
        ~isfinite(prefCohNoisePC) || ~isfinite(probeCohNoisePC) || probeCohNoisePC <= 0
    error('updateAcrossOffsetSummaries:BadNoiseAmplitude', ...
        'Invalid pref/probe coherence noise amplitudes: pref=%g, probe=%g.', ...
        prefCohNoisePC, probeCohNoisePC);
end

nYokedProbeStreams = probeStreamCountFromSessionProbeHeader(sessionProbeHeader);
combinedProbeCohNoisePC = nYokedProbeStreams * probeCohNoisePC;

normInfo = struct();
normInfo.prefCohNoisePC = prefCohNoisePC;
normInfo.probeCohNoisePC = probeCohNoisePC;
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

probeDirDeg = abs(double(localExtractValue(sessionProbeHeader.probeDirDeg)));

if probeDirDeg > 0 && probeDirDeg < 180
    n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
    n = 1;
else
    error('updateAcrossOffsetSummaries:UnsupportedProbeDir', ...
        'Unsupported probeDirDeg for probe normalization: %g.', probeDirDeg);
end
end