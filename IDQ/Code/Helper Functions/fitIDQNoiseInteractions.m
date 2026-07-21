function interactionFits = fitIDQNoiseInteractions(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDQNoiseInteractions
%
% Fit symmetry-reduced interactions among the three change-side IDQ noise
% streams. The input predictors retain their physical coherence signs.
%
% Required trialTable variables:
%   sessionIndex, stepCoh, correct, hasStepNoise, dirIndex,
%   noisePredDir1, noisePredDir2, noisePredDir3
%
% Models:
%   main: cEff = cStep + gD*nD + gN*(nPlus+nMinus)
%   driftNonDriftInteraction: main + hDN*nD*(nPlus+nMinus)
%   nonDriftInteraction:      main + hNN*nPlus*nMinus
%   bothInteractions:         main + hDN term + hNN term
%
% Product predictors are mean-centered within session and divided by the
% pooled geometric mean of the component predictor SDs. This leaves them
% in coherence units and makes their fitted coefficients gain-like.

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "noisePredDir1", "noisePredDir2", "noisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQNoiseInteractions:MissingVariables', ...
    'trialTable is missing required variables: %s', strjoin(missingVars, ', '));
end

idx = trialTable.hasStepNoise;
sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
dirIndex = double(trialTable.dirIndex(idx));
Xabs = double([trialTable.noisePredDir1(idx), ...
  trialTable.noisePredDir2(idx), trialTable.noisePredDir3(idx)]);

valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & ...
  isfinite(dirIndex) & all(isfinite(Xabs), 2) & ...
  dirIndex >= 1 & dirIndex <= 3 & dirIndex == round(dirIndex);

sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
dirIndex = dirIndex(valid);
Xabs = Xabs(valid, :);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;

sessionThreshold = nan(size(sessionIndex));
for iSession = 1:height(sessionFits)
  idxSession = sessionIndex == sessionFits.sessionIndex(iSession);
  sessionThreshold(idxSession) = sessionFits.threshold(iSession);
end
if any(~isfinite(sessionThreshold))
  error('fitIDQNoiseInteractions:MissingSessionThreshold', ...
    'Could not assign threshold to all noisy trials.');
end

sessionAlpha = idqWeibullAlphaForThreshold(targetPerformance, ...
  sessionThreshold, betaWeibull, lapse);

nTrials = numel(correct);
row = (1:nTrials)';
plusDir = mod(dirIndex, 3) + 1;
minusDir = mod(dirIndex - 2, 3) + 1;

nD = Xabs(sub2ind(size(Xabs), row, dirIndex));
nPlus = Xabs(sub2ind(size(Xabs), row, plusDir));
nMinus = Xabs(sub2ind(size(Xabs), row, minusDir));
nN = nPlus + nMinus;

[iDN, scaleDN, rawMeanDN] = makeInteractionPredictor(nD, nN, sessionIndex);
[iNN, scaleNN, rawMeanNN] = makeInteractionPredictor(nPlus, nMinus, sessionIndex);

common = struct();
common.stepCoh = stepCoh;
common.correct = correct;
common.sessionAlpha = sessionAlpha;
common.betaWeibull = betaWeibull;
common.lapse = lapse;

Xmain = [nD, nN];
main = fitGainModel('main', ["g_D", "g_N"], Xmain, [1; 0], common);
fitDN = fitGainModel('driftNonDriftInteraction', ...
  ["g_D", "g_N", "h_DN"], [Xmain, iDN], [main.gain; 0], common);
fitNN = fitGainModel('nonDriftInteraction', ...
  ["g_D", "g_N", "h_NN"], [Xmain, iNN], [main.gain; 0], common);
fitBoth = fitGainModel('bothInteractions', ...
  ["g_D", "g_N", "h_DN", "h_NN"], [Xmain, iDN, iNN], ...
  [main.gain; 0; 0], common);

fitDN = addNestedComparison(fitDN, main);
fitNN = addNestedComparison(fitNN, main);
fitBoth = addNestedComparison(fitBoth, main);

interactionFits = struct();
interactionFits.main = main;
interactionFits.driftNonDriftInteraction = fitDN;
interactionFits.nonDriftInteraction = fitNN;
interactionFits.bothInteractions = fitBoth;
interactionFits.stratifiedMain = fitStratifiedMain(nD, nN, common);

interactionFits.scaling = struct();
interactionFits.scaling.sessionIndex = unique(sessionIndex);
interactionFits.scaling.driftNonDrift = scaleDN;
interactionFits.scaling.nonDriftPair = scaleNN;
interactionFits.scaling.rawMeanDriftNonDriftProductBySession = rawMeanDN;
interactionFits.scaling.rawMeanNonDriftProductBySession = rawMeanNN;
interactionFits.scaling.definition = ...
  ['Each raw product is mean-centered within session and divided by ' ...
   'sqrt(pooled SD(x1)*pooled SD(x2)); ' ...
   'the resulting interaction predictor has coherence units.'];

interactionFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nD, nPlus, nMinus, nN, iDN, iNN, main.pPred, fitBoth.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nD', 'nPlus', 'nMinus', 'nN', 'iDN', 'iNN', 'pMain', 'pBoth'});
interactionFits.nTrials = nTrials;
interactionFits.betaWeibull = betaWeibull;
interactionFits.lapse = lapse;
interactionFits.targetPerformance = targetPerformance;

end

%% ------------------------------------------------------------------------
function [xInteraction, scale, rawMeanBySession] = makeInteractionPredictor(x1, x2, sessionIndex)

sd1 = std(x1, 'omitnan');
sd2 = std(x2, 'omitnan');
scale = sqrt(sd1 * sd2);
if ~isfinite(scale) || scale <= 0
  error('fitIDQNoiseInteractions:BadInteractionScale', ...
    'Cannot scale an interaction predictor with zero or nonfinite SD.');
end

rawProduct = x1 .* x2;
centeredProduct = rawProduct;
sessions = unique(sessionIndex);
rawMeanBySession = nan(numel(sessions), 1);
for iSession = 1:numel(sessions)
  idx = sessionIndex == sessions(iSession);
  rawMeanBySession(iSession) = mean(rawProduct(idx), 'omitnan');
  centeredProduct(idx) = rawProduct(idx) - rawMeanBySession(iSession);
end
xInteraction = centeredProduct ./ scale;

end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName, predictorNames, X, gain0, common)

nParameters = size(X, 2);
lb = -5 * ones(nParameters, 1);
ub =  5 * ones(nParameters, 1);
gain0 = gain0(:);

objective = @(gain) negLogLikelihoodGain(gain, X, common.stepCoh, ...
  common.correct, common.sessionAlpha, common.betaWeibull, common.lapse);

opts = optimoptions('fmincon', 'Display', 'off', ...
  'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-8, ...
  'StepTolerance', 1e-10, 'MaxFunctionEvaluations', 5000, ...
  'MaxIterations', 2000);

[gainHat, nll, exitflag] = fmincon(objective, gain0, [], [], [], [], ...
  lb, ub, [], opts);
gainHat = gainHat(:);

hessian = finiteDifferenceHessian(objective, gainHat);
SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);
covariance = nan(nParameters);
if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;
  eigenvalues = eig(hessian);
  if all(eigenvalues > 0)
    covariance = inv(hessian);
    variance = diag(covariance);
    validVariance = isfinite(variance) & variance > 0;
    SE(validVariance) = sqrt(variance(validVariance));
    CI95(validVariance, :) = gainHat(validVariance) + ...
      [-1 1] .* (1.96 * SE(validVariance));
  end
end

effectiveCoh = max(common.stepCoh + X * gainHat, 0);
pPred = idqWeibullP(effectiveCoh, common.sessionAlpha, ...
  common.betaWeibull, common.lapse);

fit = struct();
fit.model = modelName;
fit.predictorNames = cellstr(predictorNames(:)');
fit.gain = gainHat;
fit.SE = SE;
fit.CI95 = CI95;
fit.covariance = covariance;
fit.hessian = hessian;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.nParameters = nParameters;
fit.nTrials = numel(common.correct);
fit.nEffectiveCohClipped = sum(common.stepCoh + X * gainHat < 0);
fit.effectiveCoh = effectiveCoh;
fit.pPred = pPred;

end

%% ------------------------------------------------------------------------
function fit = addNestedComparison(fit, reducedFit)

fit.deltaNLL = reducedFit.negLogLikelihood - fit.negLogLikelihood;
fit.likelihoodRatio = 2 * fit.deltaNLL;
fit.deltaParameters = fit.nParameters - reducedFit.nParameters;
fit.pValueLRT = gammainc(max(fit.likelihoodRatio, 0) / 2, ...
  fit.deltaParameters / 2, 'upper');

end

%% ------------------------------------------------------------------------
function stratified = fitStratifiedMain(nD, nN, common)

edges = quantile(nD, [0, 1/3, 2/3, 1]);
edges(1) = -Inf;
edges(end) = Inf;

nStrata = 3;
stratum = (1:nStrata)';
label = ["Low n_D"; "Middle n_D"; "High n_D"];
nTrials = nan(nStrata, 1);
medianND = nan(nStrata, 1);
gD = nan(nStrata, 1);
gN = nan(nStrata, 1);
gNSE = nan(nStrata, 1);
gNCI95Low = nan(nStrata, 1);
gNCI95High = nan(nStrata, 1);
NLL = nan(nStrata, 1);

for iStratum = 1:nStrata
  if iStratum < nStrata
    idx = nD >= edges(iStratum) & nD < edges(iStratum + 1);
  else
    idx = nD >= edges(iStratum) & nD <= edges(iStratum + 1);
  end

  commonThis = subsetCommon(common, idx);
  fit = fitGainModel(sprintf('mainStratum%d', iStratum), ...
    ["g_D", "g_N"], [nD(idx), nN(idx)], [1; 0], commonThis);

  nTrials(iStratum) = sum(idx);
  medianND(iStratum) = median(nD(idx), 'omitnan');
  gD(iStratum) = fit.gain(1);
  gN(iStratum) = fit.gain(2);
  gNSE(iStratum) = fit.SE(2);
  gNCI95Low(iStratum) = fit.CI95(2, 1);
  gNCI95High(iStratum) = fit.CI95(2, 2);
  NLL(iStratum) = fit.negLogLikelihood;
end

stratified = table(stratum, label, nTrials, medianND, gD, gN, gNSE, ...
  gNCI95Low, gNCI95High, NLL);

end

%% ------------------------------------------------------------------------
function commonOut = subsetCommon(common, idx)

commonOut = common;
commonOut.stepCoh = common.stepCoh(idx);
commonOut.correct = common.correct(idx);
commonOut.sessionAlpha = common.sessionAlpha(idx);

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihoodGain(gain, X, stepCoh, correct, ...
  sessionAlpha, betaWeibull, lapse)

effectiveCoh = max(stepCoh + X * gain(:), 0);
p = idqWeibullP(effectiveCoh, sessionAlpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(objective, x)

x = x(:);
nParameters = numel(x);
H = nan(nParameters);
h = 3e-4 .* max(1, abs(x));
f0 = objective(x);

for i = 1:nParameters
  ei = zeros(nParameters, 1);
  ei(i) = h(i);
  fp = objective(x + ei);
  fm = objective(x - ei);
  H(i, i) = (fp - 2*f0 + fm) / h(i)^2;

  for j = i+1:nParameters
    ej = zeros(nParameters, 1);
    ej(j) = h(j);
    fpp = objective(x + ei + ej);
    fpm = objective(x + ei - ej);
    fmp = objective(x - ei + ej);
    fmm = objective(x - ei - ej);
    H(i, j) = (fpp - fpm - fmp + fmm) / ...
      (4 * h(i) * h(j));
    H(j, i) = H(i, j);
  end
end

end
