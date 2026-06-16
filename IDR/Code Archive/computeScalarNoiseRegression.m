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
%   chosenSidesAll        nTrials vector; 1-based changed patch index
%   changeSidesAll        nTrials vector; 1-based changed patch index
%   changeIndicesAll      nTrials vector; 1 = DEC, 2 = INC
%
% Sign convention
%   DchangePref   =  stepSign * change-side preferred noise
%   DnoChangePref = -stepSign * no-change-side preferred noise
%   DchosenPref   =  stepSign * chosen-side preferred noise
%   DnotChosenPref = -stepSign * not-chosen-side preferred noise
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

load(probeSessionPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', ...
    'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', ...
    'chosenSidesAll', 'changeIndicesAll');

% The current primary analysis is 0-250 ms post-step. In these matrices the
% last nStepFrames are the step period. A nonzero WinStartMS is not yet
% supported because it would require specifying a subwindow within the step.
assert(winStartMS == 0, 'Only WinStartMS == 0 is currently supported.');
assert(winDurMS <= sessionHeader.stepMS, 'WinDurMS cannot exceed sessionHeader.stepMS.');

meanPrefByPatch = scalarNoiseWindowMean(prefNoiseByPatch, sessionHeader, winStartMS, winDurMS);

[nPatches, nTrials] = size(meanPrefByPatch);
assert(nPatches == 2, 'Expected two patches after averaging preferred noise.');

correct = trialOutcomesAll(:) == 0;
chosenSidesAll = chosenSidesAll(:);
changeSidesAll = changeSidesAll(:);
changeIndicesAll = changeIndicesAll(:);

assert(numel(correct) == nTrials, 'trialOutcomesAll length does not match noise trials.');
assert(numel(chosenSidesAll) == nTrials, 'chosenSidesAll length does not match noise trials.');
assert(numel(changeSidesAll) == nTrials, 'changeSidesAll length does not match noise trials.');
assert(numel(changeIndicesAll) == nTrials, 'changeIndicesAll length does not match noise trials.');

[noiseChangePref, noiseNoChangePref] = changeNoChangeNoise(meanPrefByPatch, changeSidesAll);
[noiseChosenPref, noiseNotChosenPref] = chosenNotChosenNoise(meanPrefByPatch, chosenSidesAll);
stepSign = scalarNoiseStepSign(changeIndicesAll);

meanProbeByPatch = scalarNoiseWindowMean(probeNoiseByPatch, sessionHeader, winStartMS, winDurMS);

[noiseChangeProbe, noiseNoChangeProbe] = changeNoChangeNoise(meanProbeByPatch, changeSidesAll);
[noiseChosenProbe, noiseNotChosenProbe] = chosenNotChosenNoise(meanProbeByPatch, chosenSidesAll);

DchangePref =  stepSign .* noiseChangePref;
DnoChangePref = -stepSign .* noiseNoChangePref;
DchosenPref =  stepSign .* noiseChosenPref;
DnotChosenPref = -stepSign .* noiseNotChosenPref;
DdiffPref = DchangePref + DnoChangePref;
DchangeProbe   =  stepSign .* noiseChangeProbe;
DnoChangeProbe = -stepSign .* noiseNoChangeProbe;
DchosenProbe   =  stepSign .* noiseChosenProbe;
DnotChosenProbe = -stepSign .* noiseNotChosenProbe;
DdiffProbe     =  DchangeProbe + DnoChangeProbe;

ZchangePref = localZ(DchangePref);
ZnoChangePref = localZ(DnoChangePref);
ZchosenPref = localZ(DchosenPref);
ZnotChosenPref = localZ(DnotChosenPref);
ZdiffPref = localZ(DdiffPref);
ZchangeProbe   = localZ(DchangeProbe);
ZnoChangeProbe = localZ(DnoChangeProbe);
ZchosenProbe   = localZ(DchosenProbe);
ZnotChosenProbe = localZ(DnotChosenProbe);
ZdiffProbe     = localZ(DdiffProbe);

valid = isfinite(correct) & ...
        isfinite(DchangePref) & isfinite(DnoChangePref) & ...
        isfinite(DchosenPref) & isfinite(DnotChosenPref) & ...
        isfinite(DdiffPref) & ...
        isfinite(DchangeProbe) & isfinite(DnoChangeProbe) & ...
        isfinite(DchosenProbe) & isfinite(DnotChosenProbe) & ...
        isfinite(DdiffProbe);

T = table();
T.correct = logical(correct(valid));
T.chosenSide = chosenSidesAll(valid);
T.changeSide = changeSidesAll(valid);
T.changeIndex = changeIndicesAll(valid);
T.isInc = changeIndicesAll(valid) == 2;
T.stepSign = stepSign(valid);

T.noiseChosenPref = noiseChosenPref(valid);
T.noiseNotChosenPref = noiseNotChosenPref(valid);
T.noiseChangePref = noiseChangePref(valid);
T.noiseNoChangePref = noiseNoChangePref(valid);

T.DchangePref = DchangePref(valid);
T.DnoChangePref = DnoChangePref(valid);
T.DchosenPref = DchosenPref(valid);
T.DnotChosenPref = DnotChosenPref(valid);
T.DdiffPref = DdiffPref(valid);

T.ZchosenPref = ZchosenPref(valid);
T.ZnotChosenPref = ZnotChosenPref(valid);
T.ZchangePref = ZchangePref(valid);
T.ZnoChangePref = ZnoChangePref(valid);
T.ZdiffPref = ZdiffPref(valid);

T.DchosenProbe = DchosenProbe(valid);
T.DnotChosenProbe = DnotChosenProbe(valid);
T.DchangeProbe = DchangeProbe(valid);
T.DnoChangeProbe = DnoChangeProbe(valid);
T.DdiffProbe = DdiffProbe(valid);

T.ZchangeProbe = ZchangeProbe(valid);
T.ZnoChangeProbe = ZnoChangeProbe(valid);
T.ZchosenProbe = ZchosenProbe(valid);
T.ZnotChosenProbe = ZnotChosenProbe(valid);
T.ZdiffProbe = ZdiffProbe(valid);

models = struct();
models.chosenPref = fitglm(T, 'cho ~ DchosenPref', 'Distribution', 'binomial', 'Link', 'logit');
models.notChosenPref = fitglm(T, 'correct ~ DnotChosenPref', 'Distribution', 'binomial', 'Link', 'logit');
models.changePref = fitglm(T, 'correct ~ DchangePref', 'Distribution', 'binomial', 'Link', 'logit');
models.noChangePref = fitglm(T, 'correct ~ DnoChangePref', 'Distribution', 'binomial', 'Link', 'logit');
models.bothPref = fitglm(T, 'correct ~ DchangePref + DnoChangePref', 'Distribution', 'binomial', 'Link', 'logit');
models.diffPref = fitglm(T, 'correct ~ DdiffPref', 'Distribution', 'binomial', 'Link', 'logit');

% Standardized models are for comparing effect magnitudes across sessions.
models.zChosenPref = fitglm(T, 'correct ~ ZchosenPref', 'Distribution', 'binomial', 'Link', 'logit');
models.zNotChosenPref = fitglm(T, 'correct ~ ZnotChosenPref', 'Distribution', 'binomial', 'Link', 'logit');
models.zChangePref = fitglm(T, 'correct ~ ZchangePref',  'Distribution', 'binomial', 'Link', 'logit');
models.zNoChangePref = fitglm(T, 'correct ~ ZnoChangePref', 'Distribution', 'binomial', 'Link', 'logit');
models.zBothPref = fitglm(T, 'correct ~ ZchangePref + ZnoChangePref', 'Distribution', 'binomial', 'Link', 'logit');
models.zDiffPref = fitglm(T, 'correct ~ ZdiffPref', 'Distribution', 'binomial', 'Link', 'logit');

% Probe-only models.
models.chosenProbe = fitglm(T, 'correct ~ DchosenProbe', 'Distribution', 'binomial', 'Link', 'logit');
models.notChosenProbe = fitglm(T, 'correct ~ DnotChosenProbe', 'Distribution', 'binomial', 'Link', 'logit');
models.changeProbe = fitglm(T, 'correct ~ DchangeProbe', 'Distribution', 'binomial', 'Link', 'logit');
models.noChangeProbe = fitglm(T, 'correct ~ DnoChangeProbe', 'Distribution', 'binomial', 'Link', 'logit');
models.bothProbe = fitglm(T, 'correct ~ DchangeProbe + DnoChangeProbe', 'Distribution', 'binomial', 'Link', 'logit');
models.diffProbe = fitglm(T, 'correct ~ DdiffProbe', 'Distribution', 'binomial', 'Link', 'logit');

% Preferred plus probe, chosen side only.
models.chosenPrefProbe = fitglm(T, 'correct ~ DchosenPref + DchosenProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');
models.zChosenPrefProbe = fitglm(T, 'correct ~ ZchosenPref + ZchosenProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');

% Preferred plus probe, change side only.
models.changePrefProbe = fitglm(T, 'correct ~ DchangePref + DchangeProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');
models.zChangePrefProbe = fitglm(T, 'correct ~ ZchangePref + ZchangeProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');

% Symmetric diff model.
models.diffPrefProbe = fitglm(T, 'correct ~ DdiffPref + DdiffProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');
models.zDiffPrefProbe = fitglm(T, 'correct ~ ZdiffPref + ZdiffProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');

% Full low-dimensional two-patch model.
models.fullPrefProbe = fitglm(T, ...
  'correct ~ DchangePref + DnoChangePref + DchangeProbe + DnoChangeProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');
models.zFullPrefProbe = fitglm(T, 'correct ~ ZchangePref + ZnoChangePref + ZchangeProbe + ZnoChangeProbe', ...
  'Distribution', 'binomial', 'Link', 'logit');

% What follows is selective for just the chosen data
bPref  = localModelBeta(models.chosenPrefProbe, 'DchosenPref');
bProbe = localModelBeta(models.chosenPrefProbe, 'DchosenProbe');
bPrefZ  = localModelBeta(models.zChosenPrefProbe, 'ZchosenPref');
bProbeZ = localModelBeta(models.zChosenPrefProbe, 'ZchosenProbe');

reg = struct();
reg.coeffRatios.chosenProbeOverPref_z = bProbeZ / bPrefZ;
reg.coeffRatios.chosenProbeOverPref_raw = bProbe / bPref;
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
reg.params.chosenIndicesCoding = '0=RF, 2=OPPOSITE';
reg.params.chosenSidesCoding = '0-based chosen patch index';

reg.trialTable = T;
reg.models = models;
reg.nTrials = height(T);
reg.meanCorrect = mean(T.correct);

reg.bins.chosenPref = equalCountBins(T.DchosenPref, T.correct, nBins);
reg.bins.notChosenPref = equalCountBins(T.DnotChosenPref, T.correct, nBins);
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
% changeSidesAll may be:
%   0/1 = zero-based patch index
%   1/2 = MATLAB one-based patch index
%
% Internally, changePatch is converted to MATLAB 1-based indexing.

changePatchRaw = changeSidesAll(:);
validValues = unique(changePatchRaw(isfinite(changePatchRaw)));

if all(ismember(validValues, [0 1]))
    changePatch = changePatchRaw + 1;   % convert 0/1 to 1/2
elseif all(ismember(validValues, [1 2]))
    changePatch = changePatchRaw;       % already MATLAB indexing
else
    error('changeSidesAll must contain either 0/1 or 1/2 patch indices.');
end

nTrials = size(meanNoiseByPatch, 2);
assert(numel(changePatch) == nTrials, ...
    'changeSidesAll length does not match nTrials.');

noChangePatch = 3 - changePatch;

noiseChange = nan(nTrials, 1);
noiseNoChange = nan(nTrials, 1);

for iTrial = 1:nTrials
    noiseChange(iTrial)   = meanNoiseByPatch(changePatch(iTrial), iTrial);
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

effectSummary.changeProbe = localOneModelEffect(models.changeProbe, 'DchangeProbe');
effectSummary.noChangeProbe = localOneModelEffect(models.noChangeProbe, 'DnoChangeProbe');
effectSummary.changePrefProbe_DchangePref = localOneModelEffect(models.changePrefProbe, 'DchangePref');
effectSummary.changePrefProbe_DchangeProbe = localOneModelEffect(models.changePrefProbe, 'DchangeProbe');
effectSummary.zChangePrefProbe_ZchangePref = localOneModelEffect(models.zChangePrefProbe, 'ZchangePref');
effectSummary.zChangePrefProbe_ZchangeProbe = localOneModelEffect(models.zChangePrefProbe, 'ZchangeProbe');

effectSummary.diffPrefProbe_DdiffPref = localOneModelEffect(models.diffPrefProbe, 'DdiffPref');
effectSummary.diffPrefProbe_DdiffProbe = localOneModelEffect(models.diffPrefProbe, 'DdiffProbe');
end

%% =========================================================================
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

%% =========================================================================
function beta = localModelBeta(glm, predictorName)

coef = glm.Coefficients;
row = strcmp(coef.Properties.RowNames, predictorName);

assert(any(row), 'Predictor %s not found in model.', predictorName);

beta = coef.Estimate(row);

end

% =========================================================================
function e = localOneModelEffect(glm, predictorName)
% localOneModelEffect  Extract coefficient, SE, p, CI, and odds ratio.
%
% glm is a GeneralizedLinearModel returned by fitglm.
% predictorName is the row name in glm.Coefficients.

coef = glm.Coefficients;

row = strcmp(coef.Properties.RowNames, predictorName);
assert(any(row), 'Predictor %s not found in GLM coefficients.', predictorName);

beta = coef.Estimate(row);
se   = coef.SE(row);
p    = coef.pValue(row);

ci95 = beta + [-1 1] .* 1.96 .* se;

e = struct();
e.beta = beta;
e.se = se;
e.ci95 = ci95;
e.p = p;
e.oddsRatio = exp(beta);
e.oddsRatioCI95 = exp(ci95);

end