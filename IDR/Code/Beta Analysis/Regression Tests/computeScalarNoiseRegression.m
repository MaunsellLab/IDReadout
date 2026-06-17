function reg = computeScalarNoiseRegression(probeSessionPath, varargin)
% computeScalarNoiseRegression  Scalar evidence logistic regression.
%
%   reg = computeScalarNoiseRegression(probeSessionPath)
%
% Computes mean preferred-direction coherence noise over the post-step
% window and fits low-dimensional logistic regressions predicting correct
% behavior.
%
% Required probe-session variables
%   sessionHeader
%   sessionProbeHeader
%   sideTypeNames
%   lr
%   prefNoiseByPatch      2 x nVideoFrames x nTrials
%   probeNoiseByPatch     2 x nVideoFrames x nTrials
%   trialOutcomesAll      nTrials vector; 0 means correct
%   changeSidesAll        nTrials vector; 1-based changed patch index
%   changeIndicesAll      nTrials vector; 1 = DEC, 2 = INC
%
% Sign convention
%   DchangePref   =  stepSign * change-side preferred noise
%   DnoChangePref = -stepSign * no-change-side preferred noise
%   DdiffPref     =  DchangePref + DnoChangePref
%
% where stepSign = +1 for INC and -1 for DEC. Positive predictors favor
% correct behavior.

p = inputParser;
p.addParameter('WinStartMS', 0, @isscalar);
p.addParameter('WinDurMS', 250, @isscalar);
p.addParameter('NBins', 10, @isscalar);
p.parse(varargin{:});

winStartMS = p.Results.WinStartMS;
winDurMS = p.Results.WinDurMS;
nBins = p.Results.NBins;

load(probeSessionPath, ...
    'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', ...
    'prefNoiseByPatch', 'probeNoiseByPatch', ...
    'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll');

% The current primary analysis is 0-250 ms post-step. In these matrices the
% last nStepFrames are the step period. A nonzero WinStartMS is not yet
% supported because it would require specifying a subwindow within the step.
assert(winStartMS == 0, 'Only WinStartMS == 0 is currently supported.');
assert(winDurMS <= sessionHeader.stepMS, 'WinDurMS cannot exceed sessionHeader.stepMS.');

meanPrefByPatch = scalarNoiseWindowMean(prefNoiseByPatch, sessionHeader, winStartMS, winDurMS);

[nPatches, nTrials] = size(meanPrefByPatch);
assert(nPatches == 2, 'Expected two patches after averaging preferred noise.');

correct = trialOutcomesAll(:) == 0;
changeSidesAll = changeSidesAll(:);
changeIndicesAll = changeIndicesAll(:);

assert(numel(correct) == nTrials, 'trialOutcomesAll length does not match noise trials.');
assert(numel(changeSidesAll) == nTrials, 'changeSidesAll length does not match noise trials.');
assert(numel(changeIndicesAll) == nTrials, 'changeIndicesAll length does not match noise trials.');

[noiseChangePref, noiseNoChangePref] = changeNoChangeNoise(meanPrefByPatch, changeSidesAll);
stepSign = scalarNoiseStepSign(changeIndicesAll);

DchangePref =  stepSign .* noiseChangePref;
DnoChangePref = -stepSign .* noiseNoChangePref;
DdiffPref = DchangePref + DnoChangePref;

ZchangePref = localZ(DchangePref);
ZnoChangePref = localZ(DnoChangePref);
ZdiffPref = localZ(DdiffPref);

valid = isfinite(correct) & isfinite(DchangePref) & isfinite(DnoChangePref) & ...
    isfinite(DdiffPref) & isfinite(ZchangePref) & isfinite(ZnoChangePref) & ...
    isfinite(ZdiffPref);

T = table();
T.correct = logical(correct(valid));
T.changeSide = changeSidesAll(valid);
T.changeIndex = changeIndicesAll(valid);
T.isInc = changeIndicesAll(valid) == 2;
T.stepSign = stepSign(valid);

T.noiseChangePref = noiseChangePref(valid);
T.noiseNoChangePref = noiseNoChangePref(valid);

T.DchangePref = DchangePref(valid);
T.DnoChangePref = DnoChangePref(valid);
T.DdiffPref = DdiffPref(valid);

T.ZchangePref = ZchangePref(valid);
T.ZnoChangePref = ZnoChangePref(valid);
T.ZdiffPref = ZdiffPref(valid);

models = struct();
models.changePref = fitglm(T, 'correct ~ DchangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.noChangePref = fitglm(T, 'correct ~ DnoChangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.bothPref = fitglm(T, 'correct ~ DchangePref + DnoChangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.diffPref = fitglm(T, 'correct ~ DdiffPref', ...
    'Distribution', 'binomial', 'Link', 'logit');

% Standardized models are for comparing effect magnitudes across sessions.
models.zChangePref = fitglm(T, 'correct ~ ZchangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.zNoChangePref = fitglm(T, 'correct ~ ZnoChangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.zBothPref = fitglm(T, 'correct ~ ZchangePref + ZnoChangePref', ...
    'Distribution', 'binomial', 'Link', 'logit');
models.zDiffPref = fitglm(T, 'correct ~ ZdiffPref', ...
    'Distribution', 'binomial', 'Link', 'logit');

reg = struct();
reg.probeSessionPath = probeSessionPath;
reg.sessionHeader = sessionHeader;
reg.sessionProbeHeader = sessionProbeHeader;
reg.sideTypeNames = sideTypeNames;
reg.lr = lr;

reg.params.analysisName = 'scalarNoiseRegression';
reg.params.noiseWindowName = sprintf('postStep%dto%dms', winStartMS, winStartMS + winDurMS);
reg.params.winStartMS = winStartMS;
reg.params.winDurMS = winDurMS;
reg.params.nBins = nBins;
reg.params.noiseUnits = 'mean percent coherence over window';
reg.params.primaryNoise = 'preferred direction';
reg.params.changeIndicesCoding = '1=DEC, 2=INC';
reg.params.changeSidesCoding = '1-based changed patch index';

reg.trialTable = T;
reg.models = models;
reg.nTrials = height(T);
reg.meanCorrect = mean(T.correct);

reg.bins.changePref = equalCountBins(T.DchangePref, T.correct, nBins);
reg.bins.noChangePref = equalCountBins(T.DnoChangePref, T.correct, nBins);
reg.bins.diffPref = equalCountBins(T.DdiffPref, T.correct, nBins);

reg.effectSummary = scalarNoiseEffectSummary(models);

end

% =========================================================================
function meanNoiseByPatch = scalarNoiseWindowMean(noiseByPatch, sessionHeader, winStartMS, winDurMS)
% scalarNoiseWindowMean  Mean coherence noise in a post-step window.
%
% noiseByPatch is patch x videoFrame x trial. The last nStepFrames entries
% are assumed to correspond to the full coherence-step period.

assert(ndims(noiseByPatch) == 3, 'noiseByPatch must be patch x videoFrame x trial.');
assert(size(noiseByPatch, 1) == 2, 'Expected dimension 1 to contain two patches.');

nStepFrames = round(sessionHeader.stepMS / 1000 * sessionHeader.frameRateHz);
nWinFrames = round(winDurMS / 1000 * sessionHeader.frameRateHz);
nStartFrames = round(winStartMS / 1000 * sessionHeader.frameRateHz);

nFrames = size(noiseByPatch, 2);
assert(nFrames >= nStepFrames, 'Noise matrix has fewer frames than the step window.');
assert(nStartFrames + nWinFrames <= nStepFrames, 'Requested window exceeds step period.');

stepFrame1 = nFrames - nStepFrames + 1;
winFrame1 = stepFrame1 + nStartFrames;
winFrame2 = winFrame1 + nWinFrames - 1;
winFrames = winFrame1:winFrame2;

meanNoiseByPatch = squeeze(mean(noiseByPatch(:, winFrames, :), 2, 'omitnan'));

if isvector(meanNoiseByPatch)
    meanNoiseByPatch = reshape(meanNoiseByPatch, 2, []);
end
end

% =========================================================================
function [noiseChange, noiseNoChange] = changeNoChangeNoise(meanNoiseByPatch, changeSidesAll)
% changeNoChangeNoise  Extract changed and unchanged patch noise.
%
% meanNoiseByPatch is 2 x nTrials.
% changeSidesAll is 1-based changed patch index: 1 or 2.

changePatch = changeSidesAll(:);
validValues = unique(changePatch(isfinite(changePatch)));
assert(all(ismember(validValues, [1 2])), ...
    'changeSidesAll must contain 1-based patch indices: 1 or 2.');

nTrials = size(meanNoiseByPatch, 2);
assert(numel(changePatch) == nTrials, 'changeSidesAll length does not match nTrials.');

noChangePatch = 3 - changePatch;
noiseChange = nan(nTrials, 1);
noiseNoChange = nan(nTrials, 1);

for iTrial = 1:nTrials
    noiseChange(iTrial) = meanNoiseByPatch(changePatch(iTrial), iTrial);
    noiseNoChange(iTrial) = meanNoiseByPatch(noChangePatch(iTrial), iTrial);
end
end

% =========================================================================
function stepSign = scalarNoiseStepSign(changeIndicesAll)
% scalarNoiseStepSign  Return +1 for INC and -1 for DEC.
%
% changeIndicesAll coding:
%   1 = DEC
%   2 = INC

idx = changeIndicesAll(:);
validValues = unique(idx(isfinite(idx)));
assert(all(ismember(validValues, [1 2])), ...
    'changeIndicesAll must contain 1=DEC and 2=INC.');

stepSign = nan(size(idx));
stepSign(idx == 1) = -1;
stepSign(idx == 2) =  1;
end

% =========================================================================
function z = localZ(x)

mu = mean(x, 'omitnan');
sd = std(x, 'omitnan');

if ~isfinite(sd) || sd == 0
    z = nan(size(x));
else
    z = (x - mu) ./ sd;
end
end

% =========================================================================
function bins = equalCountBins(x, y, nBins)

x = x(:);
y = double(y(:));

valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

[xx, order] = sort(x);
yy = y(order);

n = numel(xx);
edges = round(linspace(1, n + 1, nBins + 1));

bins = table();
bins.xMean = nan(nBins, 1);
bins.xSem = nan(nBins, 1);
bins.pCorrect = nan(nBins, 1);
bins.pSem = nan(nBins, 1);
bins.n = nan(nBins, 1);

for iBin = 1:nBins
    idx = edges(iBin):(edges(iBin + 1) - 1);
    if isempty(idx)
        continue;
    end

    xb = xx(idx);
    yb = yy(idx);

    bins.xMean(iBin) = mean(xb, 'omitnan');
    bins.xSem(iBin) = std(xb, 'omitnan') ./ sqrt(numel(xb));
    bins.pCorrect(iBin) = mean(yb, 'omitnan');
    bins.pSem(iBin) = sqrt(bins.pCorrect(iBin) .* (1 - bins.pCorrect(iBin)) ./ numel(yb));
    bins.n(iBin) = numel(yb);
end
end

% =========================================================================
function effectSummary = scalarNoiseEffectSummary(models)

effectSummary = struct();

effectSummary.changePref.DchangePref = oneEffect(models.changePref, 'DchangePref');
effectSummary.noChangePref.DnoChangePref = oneEffect(models.noChangePref, 'DnoChangePref');
effectSummary.diffPref.DdiffPref = oneEffect(models.diffPref, 'DdiffPref');

effectSummary.bothPref.DchangePref = oneEffect(models.bothPref, 'DchangePref');
effectSummary.bothPref.DnoChangePref = oneEffect(models.bothPref, 'DnoChangePref');

effectSummary.zBothPref.ZchangePref = oneEffect(models.zBothPref, 'ZchangePref');
effectSummary.zBothPref.ZnoChangePref = oneEffect(models.zBothPref, 'ZnoChangePref');
effectSummary.zDiffPref.ZdiffPref = oneEffect(models.zDiffPref, 'ZdiffPref');

end

% =========================================================================
function e = oneEffect(glm, predictorName)

coef = glm.Coefficients;
row = strcmp(coef.Properties.RowNames, predictorName);
assert(any(row), 'Predictor not found in GLM coefficients: %s', predictorName);

beta = coef.Estimate(row);
se = coef.SE(row);
p = coef.pValue(row);
ci95 = beta + [-1 1] .* 1.96 .* se;

e = struct();
e.beta = beta;
e.se = se;
e.ci95 = ci95;
e.p = p;
e.oddsRatio = exp(beta);
e.oddsRatioCI95 = exp(ci95);
end
