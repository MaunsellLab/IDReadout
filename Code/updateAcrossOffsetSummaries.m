function acrossOffsetSummary = updateAcrossOffsetSummaries(summaryDir, varargin)
% updateAcrossOffsetSummaries
%  MT readout fit:
%       fits a parameterized DOG readout over MT preferred direction,
%       and uses a fixed MT forward model to map that readout onto
%       predicted normalized psychophysical scale.
%
%   The fitted object is the readout a(phi), not an "effective weighting"
%   defined in kernel space.
%
% This function is intentionally modular: it does NOT construct per-session
% summaries. It loads previously saved per-session summary files, groups them
% by probe offset, applies exclusion logic, performs across-offset bootstrap
% resampling, fits the readout model, and saves a single across-offset
% summary structure.
%
%   BaselineMode handling:
%       BaselineMode = 'auto' uses a fixed readout offset b when fewer than
%       MinOffsetsForFitBaseline non-anchor offsets are available;
%       otherwise b is fit as a free parameter.
%
% REQUIRED INPUT
%   summaryDir : directory containing per-session summary .mat files
%
% NAME-VALUE OPTIONS
%   'SaveFile'        : full path to output .mat file
%   'PlotDir'         : directory for output plots
%   'NBoot'           : number of hierarchical bootstrap replicates (default 1000)
%   'CILevels'        : e.g. [68 95]
%   'Model'           : active effective-weighting model (default 'dog')
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
% -------------------------------------------------------------------------
% diag = diagnoseReadoutDOGReachability( ...
%     'TargetS45', 0.206, ...
%     'TargetS180', 0.584, ...
%     'MakePlot', true);
% 
% diag.bestParams
% diag.bestPrediction
% diag.bestDistance

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

% Legacy symbolic DOGFixedBaseline handling is not applied to the active
% readout-model fit. In the readout model, a fixed baseline must be an
% explicit numeric value for readout offset b, not a symbolic mapping from
% an empirical scale measurement such as 180 deg.

bootstrap = runHierarchicalBootstrap(offsetData, opts);
acrossOffsetSummary.bootstrap = bootstrap;

% Legacy descriptive model-fit summaries retained temporarily for backward
% compatibility, but no longer the primary interpretive fit.
modelFits = summarizeModelFits(bootstrap, empirical, opts);
acrossOffsetSummary.modelFits = modelFits;
acrossOffsetSummary.modelFitsNote = ...
    ['Legacy offset-space model summaries retained for backward ' ...
     'compatibility only. Primary interpretation should use ' ...
     'acrossOffsetSummary.readoutModel.'];

%% ---- Effective MT-to-choice weighting fit ----

offsetsDegAll   = [empirical.probeOffsetDeg];
pooledScaleAll  = [empirical.pooledScale];
bootstrapVarAll = bootstrap.fitBootstrap.offsetFitVar;

assert(numel(bootstrapVarAll) == numel(offsetsDegAll), ...
    'bootstrap offset variances do not match empirical offsets');

% Fixed MT forward model used to map the DOG readout onto predicted scale.
mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', 37.5, ...
    'phiDeg', -180:1:179);

% Store measured pooled scales by offset for downstream plotting/fits
acrossOffsetSummary.measurements = struct();
acrossOffsetSummary.measurements.offsetsDeg   = offsetsDegAll(:)';
acrossOffsetSummary.measurements.pooledScale  = pooledScaleAll(:)';
acrossOffsetSummary.measurements.bootstrapVar = bootstrapVarAll(:)';

% Fit only non-anchor offsets
isAnchor = abs(offsetsDegAll) < 1e-9;
fitOffsetsDeg = offsetsDegAll(~isAnchor);
fitScales     = pooledScaleAll(~isAnchor);
fitVars       = bootstrapVarAll(~isAnchor);

% Resolve baseline mode

nNonAnchor = numel(fitOffsetsDeg);
baselineModeRequested = lower(char(string(opts.BaselineMode)));

acrossOffsetSummary.readoutModel.activeModelName = 'dog_readout';
acrossOffsetSummary.readoutModel.sourceScaleSideType = opts.ScaleSideType;
acrossOffsetSummary.readoutModel.sourceScaleStepType = opts.ScaleStepType;
acrossOffsetSummary.readoutModel.sourceScaleMode     = opts.ScaleMode;
acrossOffsetSummary.readoutModel.mtForwardModelParams = struct( ...
    'sigmaMTDeg', mtModel.sigmaMTDeg, ...
    'phiDeg', mtModel.phiDeg);
acrossOffsetSummary.readoutModel.note = '';
acrossOffsetSummary.readoutModel = struct();
acrossOffsetSummary.readoutModel.activeModelName = 'dog_readout';
acrossOffsetSummary.readoutModel.sourceScaleSideType = opts.ScaleSideType;
acrossOffsetSummary.readoutModel.sourceScaleStepType = opts.ScaleStepType;
acrossOffsetSummary.readoutModel.sourceScaleMode     = opts.ScaleMode;
acrossOffsetSummary.readoutModel.mtForwardModelParams = struct( ...
    'sigmaMTDeg', mtModel.sigmaMTDeg, ...
    'phiDeg', mtModel.phiDeg);

% Fit only non-anchor offsets
isAnchor = abs(offsetsDegAll) < 1e-9;
fitOffsetsDeg = offsetsDegAll(~isAnchor);
fitScales     = pooledScaleAll(~isAnchor);
fitVars       = bootstrapVarAll(~isAnchor);

% Resolve BaselineMode for readout offset b only
nNonAnchor = numel(fitOffsetsDeg);
baselineModeRequested = lower(char(string(opts.BaselineMode)));

switch baselineModeRequested
    case 'auto'
        if nNonAnchor < opts.MinOffsetsForFitBaseline
            baselineModeResolved = 'fixed';
        else
            baselineModeResolved = 'fit';
        end

    case 'fixed'
        baselineModeResolved = 'fixed';

    case 'fit'
        baselineModeResolved = 'fit';

    otherwise
        error('Unknown BaselineMode: %s', char(string(opts.BaselineMode)));
end

% Resolve fixed readout baseline b if needed
fixedReadoutBaseline = [];
invalidSymbolicFixedBaseline = false;

if strcmpi(baselineModeResolved, 'fixed')
    fixedReadoutBaseline = opts.FixedReadoutBaseline;

    if ischar(fixedReadoutBaseline) || isstring(fixedReadoutBaseline)
        invalidSymbolicFixedBaseline = true;
        fixedReadoutBaseline = [];
    end
end

% With too few non-anchor offsets, the readout DOG is underconstrained.
% In that regime we only allow fitting if the user has supplied a numeric
% fixed baseline for b.
insufficientOffsetsForStableReadoutFit = ...
    nNonAnchor < opts.MinOffsetsForFitBaseline;

hasNumericFixedBaseline = ...
    isnumeric(opts.FixedReadoutBaseline) && isscalar(opts.FixedReadoutBaseline) && ...
    isfinite(opts.FixedReadoutBaseline);

allowSparseFixedBaselineFit = ...
    strcmpi(baselineModeResolved, 'fixed') && hasNumericFixedBaseline;

acrossOffsetSummary.readoutModel.baselineModeRequested = baselineModeRequested;
acrossOffsetSummary.readoutModel.baselineModeResolved  = baselineModeResolved;
acrossOffsetSummary.readoutModel.fixedBaseline         = fixedReadoutBaseline;
acrossOffsetSummary.readoutModel.minOffsetsForFitBaseline = ...
    opts.MinOffsetsForFitBaseline;

acrossOffsetSummary.readoutModel.measurements = struct();
acrossOffsetSummary.readoutModel.measurements.offsetsDeg  = offsetsDegAll(:)';
acrossOffsetSummary.readoutModel.measurements.pooledScale = pooledScaleAll(:)';
acrossOffsetSummary.readoutModel.measurements.bootstrapVar = bootstrapVarAll(:)';

if numel(fitOffsetsDeg) >= 1 && ...
        all(isfinite(fitScales)) && ...
        all(isfinite(fitVars)) && ...
        all(fitVars > 0)
  
   readoutFit = fitReadoutDOGToScales( ...
    fitOffsetsDeg, fitScales, fitVars, mtModel, ...
    'Bounds', opts.Bounds);

    acrossOffsetSummary.readoutModel.note = ...
      ['Fit of a DOG readout over MT preferred direction. ' ...
     'Predicted normalized psychophysical scale is computed through the ' ...
     'fixed MT forward model. The readout offset b is not fit because it ' ...
     'is not identifiable under the mean-subtracted forward model.'];
    acrossOffsetSummary.readoutModel.fit = readoutFit;
    acrossOffsetSummary.readoutModel.nFreeParams = readoutFit.nFreeParams;
    acrossOffsetSummary.readoutModel.phiDeg = mtModel.phiDeg;
    acrossOffsetSummary.readoutModel.readoutPhiRaw = readoutFit.readoutPhiRaw;
    acrossOffsetSummary.readoutModel.readoutPhi = readoutFit.readoutPhi;
    acrossOffsetSummary.readoutModel.readoutNormalization = 'a(0)=1';
    acrossOffsetSummary.readoutModel.paramNames = readoutFit.paramNames;
    acrossOffsetSummary.readoutModel.params = readoutFit.params;
    acrossOffsetSummary.readoutModel.paramStruct = readoutFit.paramStruct;
    acrossOffsetSummary.readoutModel.predictedAtMeasuredOffsets = ...
        predictNormalizedScaleFromReadout( ...
            offsetsDegAll, mtModel, readoutFit.params);
    
    acrossOffsetSummary.readoutModel.plotOffsetsDeg = 0:1:180;
    acrossOffsetSummary.readoutModel.plotPredictedScale = ...
        predictNormalizedScaleFromReadout( ...
            acrossOffsetSummary.readoutModel.plotOffsetsDeg, ...
            mtModel, readoutFit.params);
else
    acrossOffsetSummary.readoutModel.fit = [];
    acrossOffsetSummary.readoutModel.note = ...
        ['Readout model not fit: insufficient valid non-anchor offsets or ' ...
         'invalid scale/variance inputs for DOG readout fitting.'];
    acrossOffsetSummary.readoutModel.nFreeParams = NaN;
    acrossOffsetSummary.readoutModel.phiDeg = mtModel.phiDeg;
    acrossOffsetSummary.readoutModel.readoutPhiRaw = nan(size(mtModel.phiDeg));
    acrossOffsetSummary.readoutModel.readoutPhi = nan(size(mtModel.phiDeg));
    acrossOffsetSummary.readoutModel.readoutNormalization = 'a(0)=1';
    acrossOffsetSummary.readoutModel.paramNames = {};
    acrossOffsetSummary.readoutModel.params = [];
    acrossOffsetSummary.readoutModel.paramStruct = struct();
    acrossOffsetSummary.readoutModel.predictedAtMeasuredOffsets = [];
    acrossOffsetSummary.readoutModel.plotOffsetsDeg = 0:1:180;
    acrossOffsetSummary.readoutModel.plotPredictedScale = [];
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

phiDeg = acrossOffsetSummary.readoutModel.phiDeg;
mtp = acrossOffsetSummary.readoutModel.mtForwardModelParams;
mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', mtp.sigmaMTDeg, ...
    'phiDeg', mtp.phiDeg);

% sum(computeMTSymmetrizedDeltaM(phiDeg, 0, mtModel))
% sum(computeMTSymmetrizedDeltaM(phiDeg, 45, mtModel))
% sum(computeMTSymmetrizedDeltaM(phiDeg, 180, mtModel))


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
addParameter(p, 'AngleGridDeg', 0:1:180, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
addParameter(p, 'Model', 'dog', @(x) ischar(x) || isstring(x));
addParameter(p, 'CandidateModels', {'dog'}, @(x) ischar(x) || isstring(x) || iscellstr(x) || (iscell(x) && all(cellfun(@(c) ischar(c) || isstring(c), x))));
addParameter(p, 'BaselineMode', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'FixedReadoutBaseline', [], ...
    @(x) isempty(x) || ((isnumeric(x) && isscalar(x) && isfinite(x)) || ischar(x) || isstring(x)));
addParameter(p, 'DOGFixedBaseline', 'empirical180', ...
    @(x) (isnumeric(x) && isscalar(x) && isfinite(x)) || ischar(x) || isstring(x));
addParameter(p, 'MinOffsetsForFitBaseline', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);

parse(p, summaryDir, varargin{:});
opts = p.Results;
opts.summaryDir = summaryDir;
opts.SaveFile = char(opts.SaveFile);
opts.PlotDir  = char(opts.PlotDir);
opts.Model    = char(opts.Model);
opts.CandidateModels = normalizeModelList(opts.CandidateModels, opts.Model);
opts.FilePattern = char(opts.FilePattern);
opts.OffsetField = char(opts.OffsetField);
opts.SessionNameField = char(opts.SessionNameField);
opts.ScaleMode = char(opts.ScaleMode);
% Active option for fixed readout baseline b.
% Prefer FixedReadoutBaseline; fall back to legacy DOGFixedBaseline only if
% the new option was not supplied.
if isempty(opts.FixedReadoutBaseline)
    opts.FixedReadoutBaseline = opts.DOGFixedBaseline;
end

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
    'candidateModels', {opts.CandidateModels}, ...
    'angleUnits', 'deg', ...
    'normalization', 'scale_anchor_at_0deg', ...
    'offsetKeysDeg', [], ...
    'readoutBaselineOption', 'inactive_not_identifiable', ...
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
        curveValues(b, :) = evaluateAcrossOffsetModel(params, angleGrid, opts.Model, opts);
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
function diag = diagnoseReadoutDOGReachability(varargin)
% Brute-force diagnostic for the DOG readout model.
%
% Maps the reachable (S45, S180) pairs over a parameter grid to determine
% whether the empirical target lies within the model family's range.
%
% Example:
%   diag = diagnoseReadoutDOGReachability( ...
%       'TargetS45', 0.206, ...
%       'TargetS180', 0.584, ...
%       'MakePlot', true);

p = inputParser;
addParameter(p, 'SigmaCenterDeg', 5:5:120, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SigmaSurroundDeg', 10:5:180, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SurroundGain', 0:0.05:2.0, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'TargetS45', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'TargetS180', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MakePlot', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opt = p.Results;

mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', opt.SigmaMTDeg, ...
    'phiDeg', opt.PhiDeg);

sC = opt.SigmaCenterDeg(:)';
sS = opt.SigmaSurroundDeg(:)';
aS = opt.SurroundGain(:)';

nTotal = numel(sC) * numel(sS) * numel(aS);

S45 = nan(nTotal, 1);
S180 = nan(nTotal, 1);
paramsMat = nan(nTotal, 3);
validMask = false(nTotal, 1);

idx = 0;
for i = 1:numel(sC)
    for j = 1:numel(sS)
        if sS(j) <= sC(i)
            continue;
        end
        for k = 1:numel(aS)
            idx = idx + 1;
            params = [sC(i), sS(j), aS(k)];
            pred = predictNormalizedScaleFromReadout([45 180], mtModel, params);

            paramsMat(idx, :) = params;
            S45(idx) = pred(1);
            S180(idx) = pred(2);
            validMask(idx) = all(isfinite(pred));
        end
    end
end

paramsMat = paramsMat(1:idx, :);
S45 = S45(1:idx);
S180 = S180(1:idx);
validMask = validMask(1:idx);

paramsMat = paramsMat(validMask, :);
S45 = S45(validMask);
S180 = S180(validMask);

diag = struct();
diag.paramNames = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};
diag.params = paramsMat;
diag.S45 = S45;
diag.S180 = S180;
diag.target = [opt.TargetS45, opt.TargetS180];

if isfinite(opt.TargetS45) && isfinite(opt.TargetS180) && ~isempty(S45)
    d2 = (S45 - opt.TargetS45).^2 + (S180 - opt.TargetS180).^2;
    [bestD2, bestIdx] = min(d2);
    diag.bestIdx = bestIdx;
    diag.bestParams = paramsMat(bestIdx, :);
    diag.bestPrediction = [S45(bestIdx), S180(bestIdx)];
    diag.bestDistance = sqrt(bestD2);
else
    diag.bestIdx = [];
    diag.bestParams = [];
    diag.bestPrediction = [];
    diag.bestDistance = NaN;
end

if opt.MakePlot && ~isempty(S45)
    fig = figure; clf; hold on;
    scatter(S45, S180, 10, 'filled');
    xlabel('Predicted S(45)');
    ylabel('Predicted S(180)');
    title('Reachable (S45, S180) pairs for DOG readout model');
    box off;

    if isfinite(opt.TargetS45) && isfinite(opt.TargetS180)
        plot(opt.TargetS45, opt.TargetS180, 'kp', 'MarkerSize', 14, ...
            'MarkerFaceColor', 'k');
    end
end
end

% ========================================================================
function [aPhi, paramStruct] = evaluateReadoutDOG(phiDeg, params)
% Evaluate DOG readout over MT preferred direction.
%
% Readout:
%   a(phi) = exp(-(phi^2)/(2*sigmaC^2)) ...
%            - As * exp(-(phi^2)/(2*sigmaS^2))
%
% Constraints are enforced by the fitter, not here.

phiDeg = phiDeg(:)';

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
function predScale = predictNormalizedScaleFromReadout(offsetsDeg, mtModel, params)
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
function predScale = predictNormalizedScaleFromExplicitReadout(offsetsDeg, mtModel, aPhi)
% Predict normalized scale from an explicit readout vector a(phi).
%
% The probe template follows the stimulus construction:
%   - delta = 0 deg: single-channel template
%   - 0 < delta < 180 deg: symmetrized template for equal probe-noise
%     contributions at +delta and -delta
%   - delta = 180 deg: single-channel template

phiDeg = mtModel.phiDeg(:)';
aPhi   = aPhi(:)';

if numel(aPhi) ~= numel( ...
    phiDeg)
    error('predictNormalizedScaleFromExplicitReadout: aPhi must match mtModel.phiDeg in length.');
end

refVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, 0, mtModel));
tol = 1e-10;

offsetsDeg = offsetsDeg(:)';
predScale = nan(size(offsetsDeg));

for i = 1:numel(offsetsDeg)
    probeVal = sum(aPhi .* computeMTSymmetrizedDeltaM(phiDeg, offsetsDeg(i), mtModel));
    if isfinite(refVal) && abs(refVal) > tol
        predScale(i) = probeVal / refVal;
    else
        predScale(i) = NaN;
    end
end
end
% ========================================================================
function fitResult = fitReadoutDOGToScales(offsetsDeg, obsScale, obsVar, mtModel, varargin)
% Fit DOG readout parameters to observed non-anchor scale values under the
% project's probe-to-pref scale convention.
%
% Weighted objective:
%   sum_i (obsScale_i - predScale_i)^2 / obsVar_i

p = inputParser;
addParameter(p, 'Bounds', struct(), @(x) isstruct(x));
parse(p, varargin{:});
opt = p.Results;

offsetsDeg = offsetsDeg(:);
obsScale   = obsScale(:);
obsVar     = obsVar(:);

valid = isfinite(offsetsDeg) & isfinite(obsScale) & isfinite(obsVar) & (obsVar > 0);
offsetsDeg = offsetsDeg(valid);
obsScale   = obsScale(valid);
obsVar     = obsVar(valid);

[p0, lb, ub] = initialGuessReadoutDOG(offsetsDeg, obsScale, opt.Bounds);
obj = @(p) readoutDOGObjective(p, offsetsDeg, obsScale, obsVar, mtModel);

params = nan(size(p0));
loss = NaN;
fitSuccess = false;

try
    problem = createOptimProblem('fmincon', ...
        'objective', obj, ...
        'x0', p0, ...
        'lb', lb, ...
        'ub', ub, ...
        'options', optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point'));
    [params, loss] = fmincon(problem);
    fitSuccess = all(isfinite(params));
catch
    try
        params = fminsearch(obj, p0, optimset('Display', 'off'));
        loss = obj(params);
        fitSuccess = all(isfinite(params));
    catch
        fitSuccess = false;
    end
end

predMeasured = nan(size(offsetsDeg));
paramStruct = struct();

if fitSuccess
    predMeasured = predictNormalizedScaleFromReadout(offsetsDeg, mtModel, params).';
    [~, paramStruct] = evaluateReadoutDOG(mtModel.phiDeg, params);
end

fitResult = struct();
fitResult.modelName = 'dog_readout';
fitResult.fitSuccess = fitSuccess;
fitResult.loss = loss;
fitResult.params = params;
fitResult.paramStruct = paramStruct;
fitResult.paramNames = readoutDOGParamNames();
fitResult.nFreeParams = numel(params);
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
function plotReadoutDOGFit(fitResult, mtModel, varargin)
% Minimal plot helper for fitted readout over MT preferred direction.

p = inputParser;
addParameter(p, 'titleStr', 'DOG readout fit', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

if isempty(fitResult) || ~isfield(fitResult, 'fitSuccess') || ~fitResult.fitSuccess
    return;
end

fig = figure; clf; hold on;
plot(mtModel.phiDeg, fitResult.readoutPhi, 'k-', 'LineWidth', 2);
xlabel('\phi (deg)');
ylabel('a(\phi)');
title(char(string(p.Results.titleStr)));
box off;
end

% ========================================================================
function plotReadoutOverlapDiagnostic(acrossOffsetSummary, varargin)
% Plot fitted readout, MT templates, and their products to visualize how
% overlap determines predicted normalized scale.
%
% The key diagnostic is that scale is set by:
%   S(delta) = <a, Delta m_delta> / <a, Delta m_0>
% not by the local value a(delta) alone.

p = inputParser;
addParameter(p, 'OffsetsDeg', [0 45 180], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'titleStr', 'Fit-vs-flat overlap diagnostic', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveFile', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opt = p.Results;

rm = acrossOffsetSummary.readoutModel;
if ~isfield(rm, 'fit') || isempty(rm.fit) || ~rm.fit.fitSuccess
    return;
end

phiDeg = rm.phiDeg(:)';
aPhi   = rm.readoutPhi(:)';   % normalized display readout, a(0)=1

mtp = rm.mtForwardModelParams;
mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', mtp.sigmaMTDeg, ...
    'phiDeg', mtp.phiDeg);
offsetsDeg = opt.OffsetsDeg(:)';
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
    warning('plotReadoutOverlapDiagnostic: OffsetsDeg does not include 0. Ratios will not be shown.');
    refOverlap = NaN;
    predScale = nan(size(overlap));
else
    refOverlap = overlap(idx0);
    predScale = overlap ./ refOverlap;
end
aPhiFlat = ones(size(aPhi));
predScaleFlat = predictNormalizedScaleFromExplicitReadout(offsetsDeg, mtModel, aPhiFlat);
predScaleFit  = predScale;


% ---- Build fit-vs-flat comparison figure ----
prodTermFit  = cell(1, nOffsets);
prodTermFlat = cell(1, nOffsets);

for i = 1:nOffsets
    prodTermFit{i}  = aPhi .* deltaM{i};
    prodTermFlat{i} = ones(size(aPhi)) .* deltaM{i};   % flat readout = 1
end

fig = figure(302); clf;
t = tiledlayout(2,1, 'TileSpacing', 'compact', 'Padding', 'compact');

co = lines(nOffsets);

% ---- Top panel: readout functions and MT templates ----
nexttile; hold on;

% Readout functions
hFitReadout  = plot(phiDeg, aPhi, 'k-', 'LineWidth', 2);
hFlatReadout = plot(phiDeg, ones(size(phiDeg)), 'k--', 'LineWidth', 1.8);

% MT templates, same colors used in the lower panel
hTemplates = gobjects(1, nOffsets);
for i = 1:nOffsets
    hTemplates(i) = plot(phiDeg, deltaM{i}, '-', ...
        'Color', co(i,:), 'LineWidth', 1.5);
end

yline(0, 'k:');
xlabel('\phi (deg)');
ylabel('Value');
title(char(string(opt.titleStr)));
legend([hFitReadout, hFlatReadout, hTemplates], ...
    [{'Fitted readout a(\phi)'}, {'Flat readout = 1'}, ...
     arrayfun(@(d) sprintf('\\Delta m(\\phi;%g^\\circ)', d), offsetsDeg, ...
     'UniformOutput', false)], ...
    'Location', 'best');
box off;

% ---- Bottom panel: overlap contribution functions ----
nexttile; hold on;

legendHandles = gobjects(0);
legendLabels = {};

for i = 1:nOffsets
    h1 = plot(phiDeg, prodTermFit{i}, '-', 'Color', co(i,:), 'LineWidth', 2);
    h2 = plot(phiDeg, prodTermFlat{i}, '--', 'Color', co(i,:), 'LineWidth', 1.5);

    flatOverlap = sum(prodTermFlat{i});

    legendHandles(end+1) = h1; %#ok<AGROW>
    legendLabels{end+1} = sprintf('%g^\\circ fit, <a,\\Deltam> = %.2f', ...
        offsetsDeg(i), overlap(i)); %#ok<AGROW>

    legendHandles(end+1) = h2; %#ok<AGROW>
    legendLabels{end+1} = sprintf('%g^\\circ flat, <a,\\Deltam> = %.2f', ...
        offsetsDeg(i), flatOverlap); %#ok<AGROW>
end
yline(0, 'k:');
xlabel('\phi (deg)');
ylabel('a(\phi)\Delta m(\phi;\delta)');
legend(legendHandles, legendLabels, 'Location', 'best');
box off;
subtitle(t, ['Upper: Readout Functions and MT Templates. ' 'Lower: Contribution Functions.']);

if ~isempty(char(string(opt.SaveFile)))
    saveas(fig, char(string(opt.SaveFile)));
end

% Also print a compact numeric summary to the command window
fprintf('\nReadout overlap diagnostic:\n');
for i = 1:nOffsets
    fprintf(['  delta = %6.1f deg:  <a,Delta m> = %9.4f   ' ...
             'positive = %9.4f   negative = %9.4f'], ...
        offsetsDeg(i), overlap(i), posPart(i), negPart(i));
    if ~isnan(predScaleFit(i))
        fprintf('   S_fit(delta) = %8.4f', predScaleFit(i));
    end
    if ~isnan(predScaleFlat(i))
        fprintf('   S_flat(delta) = %8.4f', predScaleFlat(i));
    else
      fprintf('   S_flat(delta) = undefined (flat readout gives zero overlap)');
    end
    fprintf('\n');
end
end

% ========================================================================
function names = readoutDOGParamNames()
names = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};
end

% ========================================================================
function [p0, lb, ub] = initialGuessReadoutDOG(offsetsDeg, obsScale, boundsStruct)

sigmaC0 = max(10, min(50, median(offsetsDeg(offsetsDeg > 0), 'omitnan')));
if isempty(sigmaC0) || ~isfinite(sigmaC0)
    sigmaC0 = 25;
end
sigmaS0 = max(sigmaC0 + 20, 90);
As0 = 0.5;

p0 = [sigmaC0, sigmaS0, As0];
lb = [1e-3, 1e-3, 0];
ub = [300, 300, 10];

if isfield(boundsStruct, 'readoutDOG')
    B = boundsStruct.readoutDOG;
    if isfield(B, 'lb'), lb = B.lb; end
    if isfield(B, 'ub'), ub = B.ub; end
end
end

% ========================================================================
function err = readoutDOGObjective(params, offsetsDeg, obsScale, obsVar, mtModel)

pred = predictNormalizedScaleFromReadout(offsetsDeg, mtModel, params);
pred = pred(:);

resid = pred - obsScale(:);
err = sum((resid.^2) ./ obsVar(:), 'omitnan');

sigmaC = params(1);
sigmaS = params(2);
As = params(3);

if ~(sigmaS > sigmaC && sigmaC > 0 && As >= 0)
    err = err + 1e12;
end

if ~isfinite(err)
    err = 1e12;
end
end

% ========================================================================
function ang = wrapTo180Local(ang)
% Map angles to [-180, 180).

ang = mod(ang + 180, 360) - 180;
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
% NOTE:
%   Gaussian-offset support is retained only as legacy comparison
%   scaffolding. The active pathway uses DOG.
switch lower(opts.Model)
    case 'gaussian_offset'
      p0 = initialGuessGaussian(x, y);
      [lb, ub] = getGaussianBounds(opts.Bounds);
      obj = @(p) sum(w .* localResidualSquared(evaluateGaussianOffset(p, x), y), 'omitnan');
    case 'dog'
      p0 = initialGuessDOG(x, y);
      [lb, ub] = getDOGBounds(opts.Bounds);
      obj = @(p) dogObjective(p, x, y, w, opts);
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
function err = dogObjective(p, x, y, w, opts)
    yHat = evaluateDOG(p, x, opts.DOGFixedBaseline).';
    err = sum(w .* (yHat - y).^2, 'omitnan');

    sigmaC = p(1);
    sigmaS = p(2);
    minRatio = 1.25;
    if sigmaS < minRatio * sigmaC
        err = err + 1e6 + 1e3 * (minRatio * sigmaC - sigmaS);
    end
    if ~isfinite(err)
        err = 1e12;
    end
end

% ========================================================================
function y = evaluateAcrossOffsetModel(params, angleDeg, modelName, opts)

switch lower(modelName)
    case 'gaussian_offset'
        y = evaluateGaussianOffset(params, angleDeg);
    case 'dog'
      y = evaluateDOG(params, angleDeg, opts.DOGFixedBaseline);
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
function y = evaluateDOG(params, angleDeg, fixedOffset)

sigmaC = params(1);
sigmaS = params(2);
alpha  = params(3);
a      = fixedOffset;

hRaw = exp(-(angleDeg(:)'.^2) ./ (2 * sigmaC.^2)) ...
     - alpha .* exp(-(angleDeg(:)'.^2) ./ (2 * sigmaS.^2));

h0 = 1 - alpha;
if ~isfinite(h0) || abs(h0) < 1e-9
    y = nan(size(angleDeg(:)'));
    return;
end

h = hRaw ./ h0;
y = a + (1 - a) .* h;
end

% ========================================================================
function p0 = initialGuessDOG(x, y)

sigmaC0 = max(10, min(60, median(x(x > 0), 'omitnan')));
if isempty(sigmaC0) || ~isfinite(sigmaC0)
    sigmaC0 = 25;
end

sigmaS0 = max(sigmaC0 + 10, 90);
alpha0  = 0.5;

p0 = [sigmaC0, sigmaS0, alpha0];
end

% ========================================================================
function [lb, ub] = getDOGBounds(boundsStruct)

lb = [1e-3, 1e-3, 0];
ub = [300, 300, 5];

if isfield(boundsStruct, 'dog')
    B = boundsStruct.dog;
    if isfield(B, 'lb'), lb = B.lb; end
    if isfield(B, 'ub'), ub = B.ub; end
end
end

% ========================================================================
function r2 = localResidualSquared(yhat, y)

yhat = yhat(:);
y = y(:);

bad = ~isfinite(yhat) | ~isfinite(y);
r2 = (yhat - y).^2;
r2(bad) = 1e12;

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

angleGrid = bootstrap.fitBootstrap.angleGridDeg;
bootScaleMat = bootstrap.bootScaleMat;
fitWeights = [];
if isfield(bootstrap, 'fitBootstrap') && isfield(bootstrap.fitBootstrap, 'fitWeights')
    fitWeights = bootstrap.fitBootstrap.fitWeights;
end

xMeasured = [empirical.probeOffsetDeg];
xFitPoint  = [0, xMeasured];
yFitPoint  = [1, [empirical.meanScale]];

modelFits = struct();
modelFits.note = ['Legacy offset-space model summaries retained for backward ' ...
                  'compatibility only. Use acrossOffsetSummary.readoutModel ' ...
                  'for the primary DOG readout fit.'];
for iModel = 1:numel(opts.CandidateModels)
    modelName = char(opts.CandidateModels{iModel});
    modelOpts = opts;
    modelOpts.Model = modelName;

    nParams = modelParamCount(modelName);
    bootParams = nan(opts.NBoot, nParams);
    bootLoss   = nan(opts.NBoot, 1);
    fitSuccess = false(opts.NBoot, 1);
    curveValues = nan(opts.NBoot, numel(angleGrid));

    for b = 1:opts.NBoot
        scaleVec = bootScaleMat(b, :);
        yFitBoot = [1, scaleVec];
        [params, loss, ok] = fitAcrossOffsetModel(xFitPoint, yFitBoot, modelOpts, fitWeights);
        if ok
            bootParams(b, :) = params;
            bootLoss(b) = loss;
            fitSuccess(b) = true;
            curveValues(b, :) = evaluateAcrossOffsetModel(params, angleGrid, modelName, opts);
        end
    end

    good = fitSuccess;
    paramsGood = bootParams(good, :);
    curvesGood = curveValues(good, :);

    [pointEstimateParams, pointEstimateLoss, pointEstimateOK] = ...
        fitAcrossOffsetModel(xFitPoint, yFitPoint, modelOpts, fitWeights);
    if ~pointEstimateOK
        pointEstimateParams = nan(1, nParams);
        pointEstimateLoss = NaN;
    end

    entry = struct( ...
        'isFit', any(good), ...
        'paramNames', {modelParamNames(modelName)}, ...
        'pointEstimate', pointEstimateParams, ...
        'ci68', ciFromMatrix(paramsGood, 68), ...
        'ci95', ciFromMatrix(paramsGood, 95), ...
        'angleGridDeg', angleGrid, ...
        'curvePointEstimate', evaluateAcrossOffsetModel(pointEstimateParams, angleGrid, modelName, opts), ...
        'curveMedianBootstrap', nanpercentile(curvesGood, 50), ...
        'curve68', curveCI(curvesGood, 68), ...
        'curve95', curveCI(curvesGood, 95), ...
        'fitMethod', 'bootstrap_refit', ...
        'objective', 'legacy_offset_space_least_squares', ...
        'bounds', opts.Bounds, ...
        'startPointRule', 'data-driven crude initializer', ...
        'lossPointEstimate', pointEstimateLoss, ...
        'fitNotes', ['Legacy offset-space fit retained for backward compatibility; ' ...
             'not the primary DOG readout fit over MT preferred direction.'], ... 
        'bootstrapParams', bootParams, ...
        'bootstrapLoss', bootLoss, ...
        'bootstrapFitSuccess', fitSuccess );

    for c = 1:numel(opts.CILevels)
        lvl = opts.CILevels(c);
        entry.(sprintf('ci%d', round(lvl))) = ciFromMatrix(paramsGood, lvl);
        entry.(sprintf('curve%d', round(lvl))) = curveCI(curvesGood, lvl);
    end

    if size(paramsGood, 2) >= 2 && size(paramsGood, 1) >= 2
        entry.bootstrapParamCov = cov(paramsGood, 'omitrows');
        entry.bootstrapParamCorr = corrcoef(paramsGood, 'Rows', 'pairwise');
    else
        entry.bootstrapParamCov = NaN;
        entry.bootstrapParamCorr = NaN;
    end

    modelFits.(modelName) = entry;
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
    'primaryFitObject', 'DOG readout over MT preferred direction', ...
    'baselineMode', acrossOffsetSummary.readoutModel.baselineModeResolved, ...
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
    isfield(acrossOffsetSummary.readoutModel, 'fit') && ...
    ~isempty(acrossOffsetSummary.readoutModel.fit) && ...
    isfield(acrossOffsetSummary.readoutModel.fit, 'fitSuccess') && ...
    acrossOffsetSummary.readoutModel.fit.fitSuccess;

% ---- Plot 1: observed scale by offset; overlay prediction if available ----
fig1 = figure(300); clf; hold on;

hObs = errorbar(offsets, obsScale, ...
    obsScale - ci68(:,1)', ci68(:,2)' - obsScale, ...
    'ko', 'LineWidth', 1.2, 'MarkerFaceColor', 'k');

if hasReadoutFit
    rm = acrossOffsetSummary.readoutModel;
    predScale = rm.predictedAtMeasuredOffsets;
    hPred = plot(offsets, predScale, 'k-', 'LineWidth', 2);

    legend([hObs, hPred], ...
        {'Observed pooled scale (68% CI)', 'Predicted from fitted readout'}, ...
        'Location', 'best');
    title('Observed vs predicted scale');
else
    legend(hObs, {'Observed pooled scale (68% CI)'}, ...
        'Location', 'best');
    title('Observed scale by probe offset');
end

xlabel('Probe offset (deg)');
ylabel('Normalized scale');
box off;

saveas(fig1, fullfile(opts.PlotDir, 'observed_vs_predicted_scale.png'));

% ---- Plot 2: fitted readout over MT preferred direction ----
if hasReadoutFit
    rm = acrossOffsetSummary.readoutModel;

    fig2 = figure(301); clf; hold on;
    plot(rm.phiDeg, rm.readoutPhi, 'k-', 'LineWidth', 2);

    xlabel('\phi (deg)');
    ylabel('a(\phi)');
    title('Fitted DOG readout over MT preferred direction');
    box off;

    saveas(fig2, fullfile(opts.PlotDir, 'fitted_readout_phi.png'));

    plotReadoutOverlapDiagnostic( ...
    acrossOffsetSummary, ...
    'OffsetsDeg', [0 45 180], ...
    'titleStr', 'Readout / MT overlap diagnostic', ...
    'SaveFile', fullfile(opts.PlotDir, 'readout_overlap_diagnostic.png'));
end

end

% ========================================================================
function n = modelParamCount(modelName)

switch lower(modelName)
    case 'gaussian_offset'
        n = 2;
    case 'dog'
        n = 3;
    otherwise
        error('Unsupported model: %s', modelName);
end

end

% ========================================================================
function names = modelParamNames(modelName)

switch lower(modelName)
    case 'gaussian_offset'
        names = {'sigmaDeg', 'offset'};
    case 'dog'
      names = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};
    otherwise
        error('Unsupported model: %s', modelName);
end

end

% ========================================================================
function models = normalizeModelList(candidateModels, activeModel)

if ischar(candidateModels) || isstring(candidateModels)
    models = {char(candidateModels)};
else
    models = cellfun(@(c) char(string(c)), candidateModels, 'UniformOutput', false);
end

models = unique(models, 'stable');
activeModel = char(string(activeModel));
if ~any(strcmpi(models, activeModel))
    models = [{activeModel}, models];
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
