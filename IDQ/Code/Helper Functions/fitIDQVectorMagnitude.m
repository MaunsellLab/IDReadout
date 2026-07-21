function vectorFits = fitIDQVectorMagnitude(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDQVectorMagnitude
%
% Fit trialwise net-motion-vector magnitude models to IDQ choices.
% Directions are separated by 120 degrees. Predictor signs remain physical
% coherence signs.
%
% Required trialTable variables:
%   sessionIndex, stepCoh, correct, hasStepNoise, dirIndex,
%   changeNoisePredDir1-3, noChangeNoisePredDir1-3
%
% Models:
%   stepOnly:
%     cEff = cStep
%
%   vectorFixed:
%     cEff = RChange - RNoChange
%
%   vectorSharedGain:
%     cEff = cStep + gVector*((RChange-cStep) - RNoChange)
%
%   vectorSideGains:
%     cEff = cStep + gChange*(RChange-cStep) - gNoChange*RNoChange
%
%   additiveReduced:
%     cEff = cStep + gD*nChangeD
%            + gChangeN*(nChangePlus+nChangeMinus)
%            + gNoChange*(nNoChangeD+nNoChangePlus+nNoChangeMinus)
%
%   additiveSigned:
%     additiveReduced + hNoChange*sum(abs(no-change noise)), with the
%     absolute-value predictor centered within session

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "changeNoisePredDir1", "changeNoisePredDir2", ...
  "changeNoisePredDir3", "noChangeNoisePredDir1", ...
  "noChangeNoisePredDir2", "noChangeNoisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQVectorMagnitude:MissingVariables', ...
    'trialTable is missing required variables: %s', strjoin(missingVars, ', '));
end

idx = trialTable.hasStepNoise;
sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
dirIndex = double(trialTable.dirIndex(idx));
XchangeAbs = double([trialTable.changeNoisePredDir1(idx), ...
  trialTable.changeNoisePredDir2(idx), trialTable.changeNoisePredDir3(idx)]);
XnoChangeAbs = double([trialTable.noChangeNoisePredDir1(idx), ...
  trialTable.noChangeNoisePredDir2(idx), trialTable.noChangeNoisePredDir3(idx)]);

valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & ...
  isfinite(dirIndex) & all(isfinite(XchangeAbs), 2) & ...
  all(isfinite(XnoChangeAbs), 2) & dirIndex >= 1 & dirIndex <= 3 & ...
  dirIndex == round(dirIndex);
sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
dirIndex = dirIndex(valid);
XchangeAbs = XchangeAbs(valid, :);
XnoChangeAbs = XnoChangeAbs(valid, :);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;
sessionThreshold = nan(size(sessionIndex));
for iSession = 1:height(sessionFits)
  idxSession = sessionIndex == sessionFits.sessionIndex(iSession);
  sessionThreshold(idxSession) = sessionFits.threshold(iSession);
end
if any(~isfinite(sessionThreshold))
  error('fitIDQVectorMagnitude:MissingSessionThreshold', ...
    'Could not assign threshold to all noisy trials.');
end
sessionAlpha = idqWeibullAlphaForThreshold(targetPerformance, ...
  sessionThreshold, betaWeibull, lapse);

nTrials = numel(correct);
row = (1:nTrials)';
plusDir = mod(dirIndex, 3) + 1;
minusDir = mod(dirIndex - 2, 3) + 1;

nChangeD = XchangeAbs(sub2ind(size(XchangeAbs), row, dirIndex));
nChangePlus = XchangeAbs(sub2ind(size(XchangeAbs), row, plusDir));
nChangeMinus = XchangeAbs(sub2ind(size(XchangeAbs), row, minusDir));
nNoChangeD = XnoChangeAbs(sub2ind(size(XnoChangeAbs), row, dirIndex));
nNoChangePlus = XnoChangeAbs(sub2ind(size(XnoChangeAbs), row, plusDir));
nNoChangeMinus = XnoChangeAbs(sub2ind(size(XnoChangeAbs), row, minusDir));

rChange = vectorMagnitude120(stepCoh + nChangeD, nChangePlus, nChangeMinus);
rNoChange = vectorMagnitude120(nNoChangeD, nNoChangePlus, nNoChangeMinus);
xChangeVector = rChange - stepCoh;
xNoChangeVector = -rNoChange;
xSharedVector = xChangeVector + xNoChangeVector;

common = struct();
common.stepCoh = stepCoh;
common.correct = correct;
common.sessionAlpha = sessionAlpha;
common.betaWeibull = betaWeibull;
common.lapse = lapse;

stepOnly = packageFixedModel('stepOnly', "Step only", stepCoh, 0, common);
vectorFixed = packageFixedModel('vectorFixed', "Physical vector", ...
  rChange - rNoChange, 0, common);
vectorSharedGain = fitGainModel('vectorSharedGain', "g_Vector", ...
  xSharedVector, 1, common);
vectorSideGains = fitGainModel('vectorSideGains', ...
  ["g_Change", "g_NoChange"], [xChangeVector, xNoChangeVector], ...
  [1; 1], common);

XadditiveReduced = [nChangeD, nChangePlus + nChangeMinus, ...
  nNoChangeD + nNoChangePlus + nNoChangeMinus];
additiveReduced = fitGainModel('additiveReduced', ...
  ["g_D", "g_ChangeN", "g_NoChange"], XadditiveReduced, ...
  [1; -0.1; -0.1], common);
aNoChange = centerWithinSession(sum(abs(XnoChangeAbs), 2), sessionIndex);
additiveSigned = fitGainModel('additiveSigned', ...
  ["g_D", "g_ChangeN", "g_NoChange", "h_NoChange"], ...
  [XadditiveReduced, aNoChange], [additiveReduced.gain; -0.3], common);
additiveSigned = addNestedComparison(additiveSigned, additiveReduced);

vectorSharedGain = addNestedComparisonToFixed(vectorSharedGain, vectorFixed, 1);
vectorSideGains = addNestedComparison(vectorSideGains, vectorSharedGain);

models = {stepOnly, vectorFixed, vectorSharedGain, vectorSideGains, ...
  additiveReduced, additiveSigned};
for iModel = 1:numel(models)
  models{iModel}.AIC = 2 * models{iModel}.nParameters + ...
    2 * models{iModel}.negLogLikelihood;
end
stepOnly = models{1};
vectorFixed = models{2};
vectorSharedGain = models{3};
vectorSideGains = models{4};
additiveReduced = models{5};
additiveSigned = models{6};
bestAIC = min(cellfun(@(x) x.AIC, models));
for iModel = 1:numel(models)
  models{iModel}.deltaAIC = models{iModel}.AIC - bestAIC;
end
stepOnly.deltaAIC = models{1}.deltaAIC;
vectorFixed.deltaAIC = models{2}.deltaAIC;
vectorSharedGain.deltaAIC = models{3}.deltaAIC;
vectorSideGains.deltaAIC = models{4}.deltaAIC;
additiveReduced.deltaAIC = models{5}.deltaAIC;
additiveSigned.deltaAIC = models{6}.deltaAIC;

vectorFits = struct();
vectorFits.stepOnly = stepOnly;
vectorFits.vectorFixed = vectorFixed;
vectorFits.vectorSharedGain = vectorSharedGain;
vectorFits.vectorSideGains = vectorSideGains;
vectorFits.additiveReduced = additiveReduced;
vectorFits.additiveSigned = additiveSigned;
vectorFits.modelOrder = {'stepOnly', 'vectorFixed', 'vectorSharedGain', ...
  'vectorSideGains', 'additiveReduced', 'additiveSigned'};
vectorFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nChangeD, nChangePlus, nChangeMinus, nNoChangeD, nNoChangePlus, ...
  nNoChangeMinus, rChange, rNoChange, xChangeVector, xNoChangeVector, ...
  xSharedVector, aNoChange, stepOnly.pPred, vectorFixed.pPred, ...
  vectorSharedGain.pPred, vectorSideGains.pPred, additiveReduced.pPred, ...
  additiveSigned.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nChangeD', 'nChangePlus', 'nChangeMinus', 'nNoChangeD', ...
  'nNoChangePlus', 'nNoChangeMinus', 'rChange', 'rNoChange', ...
  'xChangeVector', 'xNoChangeVector', 'xSharedVector', 'aNoChange', 'pStepOnly', ...
  'pVectorFixed', 'pVectorShared', 'pVectorSide', 'pAdditiveReduced', ...
  'pAdditiveSigned'});
vectorFits.nTrials = nTrials;
vectorFits.betaWeibull = betaWeibull;
vectorFits.lapse = lapse;
vectorFits.targetPerformance = targetPerformance;
vectorFits.vectorDefinition = ...
  'R^2 = a^2+b^2+c^2-a*b-a*c-b*c for directions separated by 120 degrees.';

end

%% ------------------------------------------------------------------------
function xCentered = centerWithinSession(x, sessionIndex)

xCentered = x;
sessions = unique(sessionIndex);
for iSession = 1:numel(sessions)
  idx = sessionIndex == sessions(iSession);
  xCentered(idx) = x(idx) - mean(x(idx), 'omitnan');
end

end

%% ------------------------------------------------------------------------
function r = vectorMagnitude120(a, b, c)

rSquared = a.^2 + b.^2 + c.^2 - a.*b - a.*c - b.*c;
r = sqrt(max(rSquared, 0));

end

%% ------------------------------------------------------------------------
function fit = packageFixedModel(modelName, predictorName, effectiveCoh, nParameters, common)

effectiveCoh = max(effectiveCoh, 0);
pPred = idqWeibullP(effectiveCoh, common.sessionAlpha, ...
  common.betaWeibull, common.lapse);
pPred = min(max(pPred, eps), 1 - eps);
nll = -sum(common.correct .* log(pPred) + ...
  (1 - common.correct) .* log(1 - pPred));

fit = struct();
fit.model = modelName;
fit.predictorNames = cellstr(string(predictorName));
fit.gain = [];
fit.SE = [];
fit.CI95 = [];
fit.covariance = [];
fit.hessian = [];
fit.negLogLikelihood = nll;
fit.exitflag = NaN;
fit.nParameters = nParameters;
fit.nTrials = numel(common.correct);
fit.nEffectiveCohClipped = sum(effectiveCoh <= 0);
fit.effectiveCoh = effectiveCoh;
fit.pPred = pPred;

end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName, predictorNames, X, gain0, common)

nParameters = size(X, 2);
lb = -5 * ones(nParameters, 1);
ub =  5 * ones(nParameters, 1);
objective = @(gain) negLogLikelihoodGain(gain, X, common.stepCoh, ...
  common.correct, common.sessionAlpha, common.betaWeibull, common.lapse);
opts = optimoptions('fmincon', 'Display', 'off', ...
  'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-8, ...
  'StepTolerance', 1e-10, 'MaxFunctionEvaluations', 5000, ...
  'MaxIterations', 2000);
[gainHat, nll, exitflag] = fmincon(objective, gain0(:), ...
  [], [], [], [], lb, ub, [], opts);
gainHat = gainHat(:);

hessian = finiteDifferenceHessian(objective, gainHat);
SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);
covariance = nan(nParameters);
if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;
  if all(eig(hessian) > 0)
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
fit.predictorNames = cellstr(string(predictorNames(:)'));
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
function fit = addNestedComparisonToFixed(fit, fixedFit, deltaParameters)

fit.deltaNLL = fixedFit.negLogLikelihood - fit.negLogLikelihood;
fit.likelihoodRatio = 2 * fit.deltaNLL;
fit.deltaParameters = deltaParameters;
fit.pValueLRT = gammainc(max(fit.likelihoodRatio, 0) / 2, ...
  deltaParameters / 2, 'upper');

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
