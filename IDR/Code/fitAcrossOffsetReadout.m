function fitSummary = fitAcrossOffsetReadout(measurements, variances, offsetsDeg, varargin)
% fitAcrossOffsetReadout  Fit shared MT/readout models to normalized measurements.
%
% fitSummary = fitAcrossOffsetReadout(measurements, variances, offsetsDeg)
%
% This function is deliberately agnostic about whether measurements came
% from kernel scales or regression-beta ratios. The caller is responsible
% for preparing one normalized scalar measurement and variance per offset.
%
% Name-value options:
%   'NSessions'          Number of contributing sessions at each offset
%   'Bounds'             Optional DOG parameter bounds structure
%   'SigmaMTDeg'         MT tuning sigma (default 37.5 deg)
%   'PhiDeg'             Shared MT/readout grid (default -180:179)
%   'SourceMeasureType'  Descriptive label, e.g. 'kernelScale' or 'betaRatio'
%   'SourceSideType'     Optional source side-type description
%   'SourceStepType'     Optional source step-type description
%   'SourceMode'         Optional source measurement-mode description
%
% The 0-deg value is an analytical normalization anchor and is not fit.

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'measurements', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'variances', @(x) isnumeric(x) && isvector(x));
addRequired(p, 'offsetsDeg', @(x) isnumeric(x) && isvector(x));
addParameter(p, 'NSessions', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
addParameter(p, 'Bounds', struct(), @isstruct);
addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
addParameter(p, 'SourceMeasureType', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SourceSideType', '', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'SourceStepType', '', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'SourceMode', '', @(x) ischar(x) || isstring(x));
parse(p, measurements, variances, offsetsDeg, varargin{:});
optsIn = p.Results;

measurements = measurements(:)';
variances = variances(:)';
offsetsDeg = offsetsDeg(:)';
assert(numel(measurements) == numel(offsetsDeg), ...
    'fitAcrossOffsetReadout:SizeMismatch', 'Measurements and offsets must have equal lengths.');
assert(numel(variances) == numel(offsetsDeg), ...
    'fitAcrossOffsetReadout:SizeMismatch', 'Variances and offsets must have equal lengths.');

if isempty(optsIn.NSessions)
    nSessions = ones(size(offsetsDeg));
else
    nSessions = optsIn.NSessions(:)';
    assert(numel(nSessions) == numel(offsetsDeg), ...
        'fitAcrossOffsetReadout:SizeMismatch', ...
        'NSessions and offsets must have equal lengths.');
end

isAnchor = abs(offsetsDeg) < 1e-9;
validFitOffset = ~isAnchor & ...
    nSessions > 0 & ...
    isfinite(offsetsDeg) & ...
    isfinite(measurements) & ...
    isfinite(variances) & ...
    variances > 0;

fitOffsetsDeg = offsetsDeg(validFitOffset);
fitMeasurements = measurements(validFitOffset);
fitVars = variances(validFitOffset);

mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', optsIn.SigmaMTDeg, ...
    'phiDeg', optsIn.PhiDeg);

% Preserve the historical option field names used by the model summaries.
opts = struct();
opts.Bounds = optsIn.Bounds;
opts.ScaleSideType = optsIn.SourceSideType;
opts.ScaleStepType = optsIn.SourceStepType;
opts.ScaleMode = optsIn.SourceMode;

readoutModels = struct();
readoutModels.signedDOG = buildReadoutDOGModelSummary( ...
    'signedDOG', 'signed', ...
    fitOffsetsDeg, fitMeasurements, fitVars, ...
    offsetsDeg, measurements, variances, mtModel, opts);
readoutModels.rectifiedDOG = buildReadoutDOGModelSummary( ...
    'rectifiedDOG', 'rectified', ...
    fitOffsetsDeg, fitMeasurements, fitVars, ...
    offsetsDeg, measurements, variances, mtModel, opts);

fitSummary = struct();
fitSummary.measurements = struct( ...
    'offsetsDeg', offsetsDeg, ...
    'values', measurements, ...
    'variances', variances, ...
    'nSessions', nSessions, ...
    'validFitOffset', validFitOffset, ...
    'sourceMeasureType', char(string(optsIn.SourceMeasureType)));
fitSummary.readoutModels = readoutModels;
fitSummary.readoutModelComparison = compareReadoutDOGModels(readoutModels);
fitSummary.readoutModel = readoutModels.signedDOG;
end

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
