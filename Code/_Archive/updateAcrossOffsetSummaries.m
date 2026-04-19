function acrossOffsetSummary = updateAcrossOffsetSummaries(summaryDir, varargin)
% updateAcrossOffsetSummaries
% This function currently supports two conceptually distinct fit layers:
%   1) descriptive across-offset scale fits in scale space
%   2) mechanistic MT-readout fits using a fixed MT tuning kernel and
%      linear readout, normalized so predicted scale at 0 deg equals 1
%
% This function is intentionally modular: it does NOT construct per-session
% summaries. It loads previously saved per-session summary files, groups them
% by probe offset, applies exclusion logic, performs across-offset bootstrap
% resampling, fits a normalized Gaussian-plus-offset model (and optionally
% other models later), and saves a single across-offset summary structure.
%
% REQUIRED INPUT
%   summaryDir : directory containing per-session summary .mat files
%
% NAME-VALUE OPTIONS
%   'SaveFile'        : full path to output .mat file
%   'PlotDir'         : directory for output plots
%   'NBoot'           : number of hierarchical bootstrap replicates (default 1000)
%   'CILevels'        : e.g. [68 95]
%   'Model'           : descriptive scale-space fit model, currently 'gaussian_offset' (default)
%   'AngleGridDeg'    : dense angle grid for fitted curves (default 0:1:180)
%   'ExcludeFcn'      : function handle, [tf, reason] = f(sessionStruct)
%   'Verbose'         : logical, default false
%   'MakePlots'       : logical, default true
%   'FilePattern'     : file pattern in summaryDir (default '*.mat')
%   'OffsetField'     : field name carrying probe offset in session summary
%   'SessionNameField': field name carrying session identifier
%   'ScaleMode'       : session scale metric to use (default 'scaleFit')
%   'ScaleSideType'   : side type key/index for scale extraction (default 'change')
%   'ScaleStepType'   : step type key/index for scale extraction (default 'inc')
%   'Bounds'          : struct of optional model parameter bounds
%   'RandomSeed'      : [] or scalar seed for reproducibility
%
% REQUIRED SESSION-LEVEL CONTRACT (minimum)
%   Each per-session summary file should contain a variable/struct with:
%       header.probeDirDeg.data     (numeric scalar probe direction / offset)
%       sessionName                 (char/string)  [or inferable from filename]
%       header                      (struct)
%       compStats                   (struct)
%       hitStats                    (struct)
%       bootstrapSource            (struct)  % resampling-ready sessionData for computeSessionKernels()
%
%   The bootstrapSource should preferably contain either:
%       sessionData  (struct accepted by computeSessionKernels)
%   or the equivalent fields directly:
%       prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll,
%       changeSidesAll, changeIndicesAll, header, lr
%
%   The session-level scale estimate is then recomputed under trial
%   resampling via computeSessionKernels(sessionData, trialIdx).
%
% OUTPUT
%   acrossOffsetSummary : top-level struct saved to disk and returned
%
% NOTES
%   Descriptive scale-space fit:
%       w(theta) = (1-b)*exp(-(theta.^2)/(2*sigma^2)) + b
%   so that w(0) == 1 by construction.
%
%   Mechanistic MT-readout fit:
%       uses a fixed MT tuning kernel and linear readout, with predicted
%       scale normalized so that scale(0) == 1.
%
%   The mechanistic readout fit is currently exploratory and should not be
%   over-interpreted while the 180 deg dataset remains sparse.
%
%   Future models (e.g. DOG) can be added without changing the top-level
%   summary struct organization.
%
% -------------------------------------------------------------------------

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

modelFits = summarizeModelFits(bootstrap, empirical, opts);
acrossOffsetSummary.modelFits = modelFits;

%% ---- Mechanistic MT-readout fit (active model: gaussian_offset) ----

offsetsDegAll   = [empirical.probeOffsetDeg];
pooledScaleAll  = [empirical.pooledScale];
bootstrapVarAll = bootstrap.fitBootstrap.offsetFitVar;
assert(numel(bootstrapVarAll) == numel(offsetsDegAll), ...
    'bootstrap offset variances do not match empirical offsets');

% Fixed MT assumptions
mt = makeMTKernel( ...
    'sigmaDeg', 37.5, ...
    'nullRatioAbs', 1/3, ...
    'phiDeg', -180:1:179);

% Store measurement layer
acrossOffsetSummary.measurements = struct();
acrossOffsetSummary.measurements.offsetsDeg   = offsetsDegAll(:)';
acrossOffsetSummary.measurements.pooledScale  = pooledScaleAll(:)';
acrossOffsetSummary.measurements.bootstrapVar = bootstrapVarAll(:)';

% Fit only non-anchor offsets
isAnchor = abs(offsetsDegAll) < 1e-9;
fitOffsetsDeg = offsetsDegAll(~isAnchor);
fitScales     = pooledScaleAll(~isAnchor);
fitVars       = bootstrapVarAll(~isAnchor);

if numel(fitOffsetsDeg) >= 1 && ...
        all(isfinite(fitScales)) && ...
        all(isfinite(fitVars)) && ...
        all(fitVars > 0)

    activeModelName = 'gaussian_offset';
    readoutFit = fitReadoutModelToScales(fitOffsetsDeg, fitScales, fitVars, mt, activeModelName);

    acrossOffsetSummary.readoutModel = struct();
    acrossOffsetSummary.readoutModel.activeModelName = activeModelName;
    acrossOffsetSummary.readoutModel.mtParams = struct( ...
        'sigmaDeg', mt.sigmaDeg, ...
        'nullRatioAbs', mt.nullRatioAbs, ...
        'phiDeg', mt.phiDeg);
    acrossOffsetSummary.readoutModel.sourceScaleSideType = opts.ScaleSideType;
    acrossOffsetSummary.readoutModel.sourceScaleStepType = opts.ScaleStepType;
    acrossOffsetSummary.readoutModel.sourceScaleMode     = opts.ScaleMode;
    acrossOffsetSummary.readoutModel.note = ...
        'Exploratory mechanistic readout fit; interpret cautiously while 180 deg data remain sparse.';

    acrossOffsetSummary.readoutModel.fit = readoutFit;
    acrossOffsetSummary.readoutModel.predictedAtMeasuredOffsets = ...
        predictScalesFromReadout(offsetsDegAll, mt, activeModelName, readoutFit.params);

    acrossOffsetSummary.readoutModel.plotOffsetsDeg = 0:1:180;
    acrossOffsetSummary.readoutModel.plotPredictedScale = ...
        predictScalesFromReadout(acrossOffsetSummary.readoutModel.plotOffsetsDeg, ...
                                 mt, activeModelName, readoutFit.params);
else
    acrossOffsetSummary.readoutModel = struct();
    acrossOffsetSummary.readoutModel.activeModelName = 'gaussian_offset';
    acrossOffsetSummary.readoutModel.mtParams = struct( ...
        'sigmaDeg', mt.sigmaDeg, ...
        'nullRatioAbs', mt.nullRatioAbs, ...
        'phiDeg', mt.phiDeg);
    acrossOffsetSummary.readoutModel.sourceScaleSideType = opts.ScaleSideType;
    acrossOffsetSummary.readoutModel.sourceScaleStepType = opts.ScaleStepType;
    acrossOffsetSummary.readoutModel.sourceScaleMode     = opts.ScaleMode;
    acrossOffsetSummary.readoutModel.fit = [];
    acrossOffsetSummary.readoutModel.note = ...
        'Mechanistic readout model not fit: insufficient valid non-anchor offsets. Exploratory use only while 180 deg data remain sparse.';
end

acrossOffsetSummary.history = updateHistory(acrossOffsetSummary, opts);

saveAcrossOffsetSummary(opts, acrossOffsetSummary);

if opts.MakePlots
    try
        makeAcrossOffsetPlots(acrossOffsetSummary, opts);
    catch ME
        warning('Plot generation failed: %s', ME.message);
    end
end

if opts.Verbose
    fprintf('updateAcrossOffsetSummaries: done. Saved summary to %s\n', opts.SaveFile);
end

if isfield(acrossOffsetSummary, 'readoutModel') && ...
        isfield(acrossOffsetSummary.readoutModel, 'fit') && ...
        ~isempty(acrossOffsetSummary.readoutModel.fit)

    plotReadoutFit(acrossOffsetSummary.readoutModel.fit, mt, ...
        'titleStr', 'Across-offset MT readout fit');
end
end

% ========================================================================
function opts = parseInputs(summaryDir, varargin)

validateattributes(summaryDir, {'char','string'}, {'nonempty'}, mfilename, 'summaryDir', 1);
summaryDir = char(summaryDir);

p = inputParser;
p.FunctionName = mfilename;

addRequired(p, 'summaryDir', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveFile', fullfile(summaryDir, 'AcrossOffsetSummaries', 'IDR_acrossOffsetSummary.mat'), @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotDir',  fullfile(summaryDir, 'AcrossOffsetSummaries', 'Plots'), @(x) ischar(x) || isstring(x));
addParameter(p, 'NBoot', 1000, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'CILevels', [68 95], @(x) isnumeric(x) && isvector(x) && all(x > 0) && all(x < 100));
addParameter(p, 'Model', 'gaussian_offset', @(x) ischar(x) || isstring(x));
addParameter(p, 'AngleGridDeg', 0:1:180, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
addParameter(p, 'ExcludeFcn', [], @(x) isempty(x) || isa(x, 'function_handle'));
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
opts.summaryDir = summaryDir;
opts.SaveFile = char(opts.SaveFile);
opts.PlotDir  = char(opts.PlotDir);
opts.Model    = char(opts.Model);
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

files = dir(fullfile(opts.summaryDir, opts.FilePattern));
sessionList = repmat(emptySessionStruct(), 0, 1);

if opts.Verbose
    fprintf('  Found %d candidate session summary files.', numel(files));
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

    [tfExclude, reasonStr] = applyExcludeFcn(sessionStruct, opts.ExcludeFcn);
    sessionStruct.isExcluded = tfExclude;
    sessionStruct.excludeReason = reasonStr;
    sessionStruct.sourceFile = thisFile;

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
sessionStruct.header = getFieldOrDefault(inStruct, 'header', struct());
if isempty(fieldnames(sessionStruct.header)) && ...
        isfield(inStruct, 'kernelFile') && ...
        (ischar(inStruct.kernelFile) || isstring(inStruct.kernelFile))
    try
        K = load(char(inStruct.kernelFile), 'header');
        if isfield(K, 'header') && isstruct(K.header)
            sessionStruct.header = K.header;
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
function [tfExclude, reasonStr] = applyExcludeFcn(sessionStruct, excludeFcn)

tfExclude = sessionStruct.isExcluded;
reasonStr = sessionStruct.excludeReason;

if isempty(excludeFcn)
    return;
end

% If the summary already excludes the file, respect that.
if tfExclude
    return;
end

try
    % Preferred compatibility path for existing excludeFile(header) logic.
    out = cell(1,2);
    [out{:}] = excludeFcn(sessionStruct.header);
    tfExclude = logical(out{1});
    if numel(out) >= 2 && ~isempty(out{2})
        reasonStr = char(string(out{2}));
    elseif tfExclude
        reasonStr = 'excluded by ExcludeFcn(header)';
    end
    return;
catch
end

try
    % Backward/forward compatibility with excludeFcn(sessionStruct).
    out = cell(1,2);
    [out{:}] = excludeFcn(sessionStruct);
    tfExclude = logical(out{1});
    if numel(out) >= 2 && ~isempty(out{2})
        reasonStr = char(string(out{2}));
    elseif tfExclude
        reasonStr = 'excluded by ExcludeFcn(sessionStruct)';
    end
catch
    try
        tfExclude = logical(excludeFcn(sessionStruct.header));
        if tfExclude
            reasonStr = 'excluded by ExcludeFcn(header)';
        end
    catch
        tfExclude = logical(excludeFcn(sessionStruct));
        if tfExclude
            reasonStr = 'excluded by ExcludeFcn(sessionStruct)';
        end
    end
end

end

% ========================================================================
function acrossOffsetSummary = initializeAcrossOffsetSummary(opts, sessionList)

acrossOffsetSummary = struct();
acrossOffsetSummary.meta = struct( ...
    'analysisDate', datetime('now'), ...
    'summaryDir', opts.summaryDir, ...
    'summaryFilesUsed', {{}}, ...
    'nCandidateFiles', numel(sessionList), ...
    'scaleMetric', opts.ScaleMode, ...
    'ciType', 'bootstrap_percentile', ...
    'ciLevels', opts.CILevels(:)', ...
    'nBoot', opts.NBoot, ...
    'bootstrapType', 'hierarchical_session_trial', ...
    'fitModelDefault', opts.Model, ...
    'angleUnits', 'deg', ...
    'normalization', 'w0_equals_1', ...
    'offsetKeysDeg', [], ...
    'notes', '' );

acrossOffsetSummary.offsetData = struct([]);
acrossOffsetSummary.empirical  = struct([]);
acrossOffsetSummary.bootstrap  = struct();
acrossOffsetSummary.modelFits  = struct();
acrossOffsetSummary.history    = struct([]);

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
    idx = find([sessionList.probeOffsetDeg] == thisOffset);
    these = sessionList(idx);

    includeMask = ~[these.isExcluded];
    included = these(includeMask);

    offsetData(k).probeOffsetDeg = thisOffset;
    offsetData(k).sessionNames = {these.sessionName};
    offsetData(k).sessionDates = [these.sessionDate];
    offsetData(k).nSessions = sum(includeMask);
    offsetData(k).nSessionsTotal = numel(these);
    offsetData(k).nTrialsBySession = [included.nTrialsTotal];
    offsetData(k).nTrialsTotal = nansum(offsetData(k).nTrialsBySession);
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

nOffsets = numel(offsetData);
nBoot    = opts.NBoot;
angleGrid = opts.AngleGridDeg(:)';

bootScaleMat = nan(nBoot, nOffsets);

% First pass: collect bootstrap pooled scales for each offset.
for b = 1:nBoot
    for k = 1:nOffsets
        bootScaleMat(b, k) = computeOffsetPooledScale(offsetData(k), opts, true);
    end
    if opts.Verbose && (b == 1 || mod(b, max(1, floor(nBoot/10))) == 0)
        fprintf('  Bootstrap %d / %d\n', b, nBoot);
    end
end

% Fit weights across offsets: inverse bootstrap variance, with 0-degree anchor
% given effectively fixed weight.
offsetFitVar = nan(1, nOffsets);
fitWeights = nan(1, nOffsets + 1);  % include 0-degree anchor
fitWeights(1) = 1e9;
for k = 1:nOffsets
    x = bootScaleMat(:, k);
    x = x(isfinite(x));
    if numel(x) >= 2
        offsetFitVar(k) = var(x, 0);
    else
        offsetFitVar(k) = NaN;
    end
end

finiteVars = offsetFitVar(isfinite(offsetFitVar) & offsetFitVar > 0);
if isempty(finiteVars)
    varFloor = 1e-6;
else
    varFloor = max(1e-6, 0.01 * median(finiteVars));
end
for k = 1:nOffsets
    if isfinite(offsetFitVar(k))
        fitWeights(k+1) = 1 ./ max(offsetFitVar(k), varFloor);
    else
        fitWeights(k+1) = 1;
    end
end

bootParams   = nan(nBoot, modelParamCount(opts.Model));
bootLoss     = nan(nBoot, 1);
fitSuccess   = false(nBoot, 1);
curveValues  = nan(nBoot, numel(angleGrid));

% Second pass: fit each bootstrap replicate using fixed empirical weights.
for b = 1:nBoot
    scaleVec = bootScaleMat(b, :);
    xFit = [0, [offsetData.probeOffsetDeg]];
    yFit = [1, scaleVec];
    [params, loss, ok] = fitAcrossOffsetModel(xFit, yFit, opts, fitWeights);
    if ok
        bootParams(b, :) = params;
        bootLoss(b) = loss;
        fitSuccess(b) = true;
        curveValues(b, :) = evaluateAcrossOffsetModel(params, angleGrid, opts.Model);
    end
end

bootstrap = struct();
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
        bootstrap.offsetBootstrap(k).(sprintf('boot%d', round(lvl))) = percentileCI(x, lvl);
    end
end

bootstrap.fitBootstrap = struct( ...
    'angleGridDeg', angleGrid, ...
    'modelName', opts.Model, ...
    'params', bootParams, ...
    'paramNames', {modelParamNames(opts.Model)}, ...
    'fitSuccess', fitSuccess, ...
    'loss', bootLoss, ...
    'curveValues', curveValues, ...
    'fitWeights', fitWeights, ...
    'offsetFitVar', offsetFitVar, ...
    'varFloor', varFloor, ...
    'medianCurve', nanpercentile(curveValues(fitSuccess, :), 50), ...
    'boot68Curve', [nanpercentile(curveValues(fitSuccess, :), 16); nanpercentile(curveValues(fitSuccess, :), 84)], ...
    'boot95Curve', [nanpercentile(curveValues(fitSuccess, :), 2.5); nanpercentile(curveValues(fitSuccess, :), 97.5)], ...
    'medianParams', nanpercentile(bootParams(fitSuccess, :), 50), ...
    'paramCI68', ciFromMatrix(bootParams(fitSuccess, :), 68), ...
    'paramCI95', ciFromMatrix(bootParams(fitSuccess, :), 95) );

for c = 1:numel(opts.CILevels)
    lvl = opts.CILevels(c);
    bootstrap.fitBootstrap.(sprintf('paramCI%d', round(lvl))) = ciFromMatrix(bootParams(fitSuccess, :), lvl);
    bootstrap.fitBootstrap.(sprintf('curveCI%d', round(lvl))) = curveCI(curveValues(fitSuccess, :), lvl);
end

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
scaleVal = selectCompStatsEntry(compStats.scale, opts.ScaleSideType, opts.ScaleStepType);

end

% ========================================================================
function scaleVal = recomputeSessionScaleBootstrap(sessionStruct, opts)
% Recompute session scale under within-session trial resampling using the
% project's computeSessionKernels(sessionData, trialIdx) pathway.

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

nTrials = size(sessionData.prefNoiseByPatch, 3);
if isempty(nTrials) || nTrials < 1
    scaleVal = sessionStruct.scalePointEstimate;
    return;
end

trialIdx = randi(nTrials, [1 nTrials]);

[~, ~, ~, ~, compStats] = computeSessionKernels(sessionData, trialIdx);
scaleVal = selectCompStatsEntry(compStats.scale, opts.ScaleSideType, opts.ScaleStepType);

end

% ========================================================================
function [kernels, kVars, hitStats, compStats] = recomputeSessionKernelStruct(sessionStruct, opts, doTrialBootstrap)
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
function pooledScale = computeOffsetPooledScale(offsetStruct, opts, doBootstrap)
% Mirror kernelAverage: recompute session kernels, pool with inverse-variance
% weighting, then compute scale from the pooled kernels.

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

sessionKernels = cell(1, numel(sessIdx));
sessionKVars   = cell(1, numel(sessIdx));
refFrameRateHz = NaN;
refNFrames = NaN;

for j = 1:numel(sessIdx)
    src = sessions(sessIdx(j));
    [kernels, kVars, ~, ~] = recomputeSessionKernelStruct(src, opts, doBootstrap);
    sessionKernels{j} = kernels;
    sessionKVars{j}   = kVars;

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

[avgKernels, ~] = poolSessionKernels(sessionKernels, sessionKVars, refNFrames);
msPerVFrame = 1000.0 / refFrameRateHz;
[scaleMat, ~, ~, ~] = kernelScaleFit(avgKernels, msPerVFrame);
pooledScale = selectCompStatsEntry(scaleMat, opts.ScaleSideType, opts.ScaleStepType);

end

% ========================================================================
function tf = hasComputeSessionKernelFields(S)

req = {'prefNoiseByPatch','probeNoiseByPatch','trialOutcomesAll', ...
       'changeSidesAll','changeIndicesAll','header'};
tf = all(isfield(S, req));

end

% ========================================================================
function [params, loss, ok] = fitAcrossOffsetModel(offsetDeg, scaleVec, opts, fitWeights)

params = nan(1, modelParamCount(opts.Model));
loss   = NaN;
ok     = false;

if nargin < 4 || isempty(fitWeights)
    fitWeights = ones(size(offsetDeg));
end

x = offsetDeg(:);
y = scaleVec(:);
w = fitWeights(:);
valid = isfinite(x) & isfinite(y) & isfinite(w) & (w > 0);
x = x(valid);
y = y(valid);
w = w(valid);

if numel(x) < modelParamCount(opts.Model)
    return;
end

switch lower(opts.Model)
    case 'gaussian_offset'
        p0 = initialGuessGaussian(x, y);
        [lb, ub] = getGaussianBounds(opts.Bounds);
        obj = @(p) sum(w .* (evaluateGaussianOffset(p, x).' - y).^2, 'omitnan');
    otherwise
        error('Unsupported model: %s', opts.Model);
end

try
    problem = createOptimProblem('fmincon', ...
        'objective', obj, ...
        'x0', p0, ...
        'lb', lb, ...
        'ub', ub, ...
        'options', optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point'));
    [params, loss] = fmincon(problem);
    ok = all(isfinite(params));
catch
    % If Optimization Toolbox is unavailable or fmincon fails, fall back.
    try
        params = fminsearch(obj, p0, optimset('Display', 'off'));
        params(1) = max(params(1), eps);
        loss = obj(params);
        ok = all(isfinite(params));
    catch
        ok = false;
    end
end

end

% ========================================================================
function y = evaluateAcrossOffsetModel(params, angleDeg, modelName)

switch lower(modelName)
    case 'gaussian_offset'
        y = evaluateGaussianOffset(params, angleDeg);
    otherwise
        error('Unsupported model: %s', modelName);
end

end

% ========================================================================
function y = evaluateGaussianOffset(params, angleDeg)

sigmaDeg = params(1);
offset   = params(2);
y = (1 - offset) .* exp(-(angleDeg(:)'.^2) ./ (2 * sigmaDeg.^2)) + offset;

end

% ========================================================================
function p0 = initialGuessGaussian(x, y)

% Crude but usually stable.
sigma0 = max(10, min(90, median(x(x > 0), 'omitnan')));
if isempty(sigma0) || ~isfinite(sigma0)
    sigma0 = 45;
end
offset0 = min(y);
if ~isfinite(offset0)
    offset0 = 0;
end
p0 = [sigma0, offset0];

end

% ========================================================================
function [lb, ub] = getGaussianBounds(boundsStruct)

lb = [1e-3, -2];
ub = [300, 1];

if isfield(boundsStruct, 'gaussian_offset')
    B = boundsStruct.gaussian_offset;
    if isfield(B, 'lb'), lb = B.lb; end
    if isfield(B, 'ub'), ub = B.ub; end
end

end

% ========================================================================
function modelFits = summarizeModelFits(bootstrap, empirical, opts)

fb = bootstrap.fitBootstrap;
good = fb.fitSuccess;
params = fb.params(good, :);
curves = fb.curveValues(good, :);

if isempty(params)
    modelFits = struct();
    modelFits.(opts.Model) = struct( ...
        'isFit', false, ...
        'paramNames', {fb.paramNames}, ...
        'pointEstimate', nan(1, modelParamCount(opts.Model)), ...
        'ci68', nan(modelParamCount(opts.Model), 2), ...
        'ci95', nan(modelParamCount(opts.Model), 2), ...
        'angleGridDeg', fb.angleGridDeg, ...
        'curvePointEstimate', nan(1, numel(fb.angleGridDeg)), ...
        'curveMedianBootstrap', nan(1, numel(fb.angleGridDeg)), ...
        'curve68', nan(2, numel(fb.angleGridDeg)), ...
        'curve95', nan(2, numel(fb.angleGridDeg)), ...
        'fitMethod', 'bootstrap_refit', ...
        'objective', 'least_squares', ...
        'bounds', opts.Bounds, ...
        'startPointRule', 'data-driven crude initializer', ...
        'lossPointEstimate', NaN, ...
        'fitNotes', 'No successful bootstrap fits');
    return;
end

xFit = [0, [empirical.probeOffsetDeg]];
yFit = [1, [empirical.meanScale]];
fitWeights = [];
if isfield(bootstrap, 'fitBootstrap') && isfield(bootstrap.fitBootstrap, 'fitWeights')
    fitWeights = bootstrap.fitBootstrap.fitWeights;
end
[pointEstimateParams, pointEstimateLoss, pointEstimateOK] = ...
    fitAcrossOffsetModel(xFit, yFit, opts, fitWeights);
if ~pointEstimateOK
    pointEstimateParams = nan(1, modelParamCount(opts.Model));
    pointEstimateLoss = NaN;
end

modelFits = struct();
modelFits.(opts.Model) = struct( ...
    'isFit', any(good), ...
    'paramNames', {fb.paramNames}, ...
    'pointEstimate', pointEstimateParams, ...
    'ci68', ciFromMatrix(params, 68), ...
    'ci95', ciFromMatrix(params, 95), ...
    'angleGridDeg', fb.angleGridDeg, ...
    'curvePointEstimate', evaluateAcrossOffsetModel(nanpercentile(params, 50), fb.angleGridDeg, opts.Model), ...
    'curveMedianBootstrap', nanpercentile(curves, 50), ...
    'curve68', curveCI(curves, 68), ...
    'curve95', curveCI(curves, 95), ...
    'fitMethod', 'bootstrap_refit', ...
    'objective', 'least_squares', ...
    'bounds', opts.Bounds, ...
    'startPointRule', 'data-driven crude initializer', ...
    'lossPointEstimate', pointEstimateLoss, ...
    'fitNotes', 'w(0)=1 normalization enforced' );

for c = 1:numel(opts.CILevels)
    lvl = opts.CILevels(c);
    modelFits.(opts.Model).(sprintf('ci%d', round(lvl))) = ciFromMatrix(params, lvl);
    modelFits.(opts.Model).(sprintf('curve%d', round(lvl))) = curveCI(curves, lvl);
end

if size(params, 2) >= 2 && size(params, 1) >= 2
    modelFits.(opts.Model).bootstrapParamCov = cov(params, 'omitrows');
    modelFits.(opts.Model).bootstrapParamCorr = corrcoef(params, 'Rows', 'pairwise');
else
    modelFits.(opts.Model).bootstrapParamCov = NaN;
    modelFits.(opts.Model).bootstrapParamCorr = NaN;
end

end

% ========================================================================
function history = updateHistory(acrossOffsetSummary, opts)

history = struct( ...
    'runDate', datetime('now'), ...
    'nOffsets', numel(acrossOffsetSummary.empirical), ...
    'offsetsDeg', [acrossOffsetSummary.empirical.probeOffsetDeg], ...
    'nSessionsByOffset', [acrossOffsetSummary.empirical.nSessions], ...
    'modelName', opts.Model, ...
    'medianParams', acrossOffsetSummary.bootstrap.fitBootstrap.medianParams, ...
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
% Minimal first-pass plots. Expand as needed.

emp = acrossOffsetSummary.empirical;
fb  = acrossOffsetSummary.bootstrap.fitBootstrap;
mf  = acrossOffsetSummary.modelFits.(opts.Model);

offsets = [emp.probeOffsetDeg];
meanScale = [emp.meanScale];
ci68 = vertcat(emp.boot68);

fig = figure(300); clf; hold on;
errorbar(offsets, meanScale, meanScale - ci68(:,1)', ci68(:,2)' - meanScale, 'ko', 'LineWidth', 1.2);
plot(fb.angleGridDeg, mf.curveMedianBootstrap, 'k-', 'LineWidth', 2);
if isfield(mf, 'curve95')
    c95 = mf.curve95;
    plot(fb.angleGridDeg, c95(1,:), 'k--');
    plot(fb.angleGridDeg, c95(2,:), 'k--');
end
xlabel('Probe offset (deg)');
ylabel('Relative scale');
title(sprintf('Across-offset fit (%s)', strrep(opts.Model, '_', '\_')));
box off;

saveas(fig, fullfile(opts.PlotDir, sprintf('acrossOffset_%s.png', opts.Model)));

fig2 = figure(301); clf;
params = fb.params(fb.fitSuccess, :);
if size(params, 2) >= 2
    scatter(params(:,1), params(:,2), 20, 'filled');
    xlabel(fb.paramNames{1});
    ylabel(fb.paramNames{2});
    title('Bootstrap parameter cloud');
    box off;
    saveas(fig2, fullfile(opts.PlotDir, sprintf('acrossOffset_%s_paramCloud.png', opts.Model)));
end

end

% ========================================================================
function n = modelParamCount(modelName)

switch lower(modelName)
    case 'gaussian_offset'
        n = 2;
    otherwise
        error('Unsupported model: %s', modelName);
end

end

% ========================================================================
function names = modelParamNames(modelName)

switch lower(modelName)
    case 'gaussian_offset'
        names = {'sigmaDeg', 'offset'};
    otherwise
        error('Unsupported model: %s', modelName);
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
    name = fileStem(sourceFile);
end

end

% ========================================================================
function out = fileStem(f)

[~, out, ~] = fileparts(f);

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
        nTrials = nansum(hs.nTrials(:));
        return;
    end
end

if isfield(inStruct, 'hitStats') && isstruct(inStruct.hitStats)
    hs = inStruct.hitStats;
    if isfield(hs, 'nTrials') && isnumeric(hs.nTrials)
        nTrials = nansum(hs.nTrials(:));
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
function offsetDeg = extractProbeOffsetDeg(inStruct, offsetField)

offsetDeg = getFieldOrDefault(inStruct, offsetField, NaN);
if isnumeric(offsetDeg) && isscalar(offsetDeg) && isfinite(offsetDeg)
    return;
end

candidateFields = {'probeOffsetDeg','probeDirDeg','probeDirectionDeg','probeOffset','probeDirection'};
for i = 1:numel(candidateFields)
    if isfield(inStruct, candidateFields{i})
        val = inStruct.(candidateFields{i});
        if isnumeric(val) && isscalar(val) && isfinite(val)
            offsetDeg = val;
            return;
        elseif isstruct(val) && isfield(val, 'data') && isnumeric(val.data) && ~isempty(val.data)
            offsetDeg = val.data(1);
            return;
        end
    end
end

% compileKernelSessionSummary stores the kernel file path, which can be used
% to recover header.probeDirDeg.data when needed.
if isfield(inStruct, 'kernelFile') && (ischar(inStruct.kernelFile) || isstring(inStruct.kernelFile))
    try
        K = load(char(inStruct.kernelFile), 'header');
        if isfield(K, 'header') && isfield(K.header, 'probeDirDeg') && ...
                isstruct(K.header.probeDirDeg) && isfield(K.header.probeDirDeg, 'data') && ...
                isnumeric(K.header.probeDirDeg.data) && ~isempty(K.header.probeDirDeg.data)
            offsetDeg = K.header.probeDirDeg.data(1);
            return;
        end
    catch
    end
end

if isfield(inStruct, 'header') && isstruct(inStruct.header)
    hdr = inStruct.header;
    headerCandidates = {'probeDirDeg','probeOffsetDeg','probeDirectionDeg','probeOffset','probeDirection'};
    for i = 1:numel(headerCandidates)
        if isfield(hdr, headerCandidates{i})
            val = hdr.(headerCandidates{i});
            if isnumeric(val) && isscalar(val) && isfinite(val)
                offsetDeg = val;
                return;
            elseif isstruct(val) && isfield(val, 'data') && isnumeric(val.data) && ~isempty(val.data)
                offsetDeg = val.data(1);
                return;
            end
        end
    end
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
function out = nanpercentile(x, p)

if isempty(x)
    out = NaN;
    return;
end
out = prctile(x, p, 1);

end

% ========================================================================
function ci = ciFromMatrix(x, level)

if isempty(x)
    ci = NaN(0,2);
    return;
end
alpha = (100 - level) / 2;
ci = [prctile(x, alpha, 1); prctile(x, 100 - alpha, 1)]';

end

% ========================================================================
function ci = curveCI(curves, level)

if isempty(curves)
    ci = NaN(2,0);
    return;
end
alpha = (100 - level) / 2;
ci = [prctile(curves, alpha, 1); prctile(curves, 100 - alpha, 1)];

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
