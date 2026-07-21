function signedFits = fitIDQNoChangeSignedNoiseAsymmetry(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDQNoChangeSignedNoiseAsymmetry
%
% Fit a symmetry-reduced signed-noise model to the three no-change-side
% streams. Predictor signs remain physical coherence signs.
%
% Null model:
%   cEff = cStep + gAll*(n1+n2+n3)
%
% Signed model:
%   cEff = cStep + gAll*(n1+n2+n3) + hAll*A
% where A = |n1|+|n2|+|n3|, centered within session.
%
% The slope for each stream is:
%   negative noise: gNegative = gAll - hAll
%   positive noise: gPositive = gAll + hAll
%
% A max-like no-change-side readout predicts hAll < 0: positive noise
% should have a more negative behavioral gain than negative noise.
%
% A constrained opponent model fixes preferred:null sensitivity at 2.5:
%   gPositive = -kPreferred
%   gNegative =  kPreferred/2.5

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "noisePredDir1", "noisePredDir2", "noisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQNoChangeSignedNoiseAsymmetry:MissingVariables', ...
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
  error('fitIDQNoChangeSignedNoiseAsymmetry:MissingSessionThreshold', ...
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

nAll = sum(Xabs, 2);
aAll = centerWithinSession(sum(abs(Xabs), 2), sessionIndex);

common = struct();
common.stepCoh = stepCoh;
common.correct = correct;
common.sessionAlpha = sessionAlpha;
common.betaWeibull = betaWeibull;
common.lapse = lapse;

main = fitGainModel('main', "g_All", nAll, -0.1, common);
signed = fitGainModel('signedAsymmetry', ["g_All", "h_All"], ...
  [nAll, aAll], [main.gain; 0], common);
signed = addNestedComparison(signed, main);
slopes = makeSignedSlopes(signed);

preferredToNullRatio = 2.5;
linearWeight = (1 - 1 / preferredToNullRatio) / 2;
absoluteWeight = (1 + 1 / preferredToNullRatio) / 2;
xOpponentConstrained = -linearWeight * nAll - absoluteWeight * aAll;
opponentConstrained = fitGainModel('opponentRatio2p5', ...
  "k_Preferred", xOpponentConstrained, max(-slopes.gain(2), 0.1), common);
opponentConstrained.preferredToNullRatio = preferredToNullRatio;
opponentConstrained.kPreferred = opponentConstrained.gain(1);
opponentConstrained.kNull = opponentConstrained.gain(1) / preferredToNullRatio;
opponentConstrained.gNegative = opponentConstrained.kNull;
opponentConstrained.gPositive = -opponentConstrained.kPreferred;
opponentConstrained.deltaNLLToFree = ...
  opponentConstrained.negLogLikelihood - signed.negLogLikelihood;
opponentConstrained.likelihoodRatioVsFree = ...
  2 * max(opponentConstrained.deltaNLLToFree, 0);
opponentConstrained.pValueLRTVsFree = gammainc( ...
  opponentConstrained.likelihoodRatioVsFree / 2, 1 / 2, 'upper');

main.AIC = 2 * main.nParameters + 2 * main.negLogLikelihood;
opponentConstrained.AIC = 2 * opponentConstrained.nParameters + ...
  2 * opponentConstrained.negLogLikelihood;
signed.AIC = 2 * signed.nParameters + 2 * signed.negLogLikelihood;
bestAIC = min([main.AIC, opponentConstrained.AIC, signed.AIC]);
main.deltaAIC = main.AIC - bestAIC;
opponentConstrained.deltaAIC = opponentConstrained.AIC - bestAIC;
signed.deltaAIC = signed.AIC - bestAIC;

signedFits = struct();
signedFits.main = main;
signedFits.signedAsymmetry = signed;
signedFits.signedSlopes = slopes;
signedFits.opponentConstrained = opponentConstrained;
signedFits.freePreferredToNullRatio = estimatePreferredToNullRatio(signed);
signedFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nD, nPlus, nMinus, nAll, aAll, xOpponentConstrained, ...
  main.pPred, signed.pPred, opponentConstrained.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nD', 'nPlus', 'nMinus', 'nAll', 'aAll', 'xOpponentConstrained', ...
  'pMain', 'pSigned', 'pOpponentConstrained'});
signedFits.nTrials = nTrials;
signedFits.betaWeibull = betaWeibull;
signedFits.lapse = lapse;
signedFits.targetPerformance = targetPerformance;
signedFits.predictorDefinition = ...
  ['A = sum(abs(no-change stream predictor)), centered within session; ' ...
   'gPositive = gAll + hAll and gNegative = gAll - hAll.'];

end

%% ------------------------------------------------------------------------
function ratio = estimatePreferredToNullRatio(fit)

g = fit.gain(1);
h = fit.gain(2);
gNegative = g - h;
gPositive = g + h;

ratio = struct();
ratio.estimate = -gPositive / gNegative;
gradient = [2 * h / gNegative^2; -2 * g / gNegative^2];
variance = gradient' * fit.covariance * gradient;
if isfinite(variance) && variance > 0
  ratio.SE = sqrt(variance);
  ratio.CI95 = ratio.estimate + [-1 1] .* 1.96 .* ratio.SE;
else
  ratio.SE = NaN;
  ratio.CI95 = [NaN NaN];
end

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
function slopes = makeSignedSlopes(fit)

Lnegative = [1 -1];
Lpositive = [1 1];
negative = makeLinearCombination(fit, Lnegative);
positive = makeLinearCombination(fit, Lpositive);

signName = ["Negative"; "Positive"];
gain = [negative.estimate; positive.estimate];
SE = [negative.SE; positive.SE];
CI95Low = [negative.CI95(1); positive.CI95(1)];
CI95High = [negative.CI95(2); positive.CI95(2)];
slopes = table(signName, gain, SE, CI95Low, CI95High);

end

%% ------------------------------------------------------------------------
function result = makeLinearCombination(fit, L)

estimate = L * fit.gain;
variance = L * fit.covariance * L';
if isfinite(variance) && variance > 0
  SE = sqrt(variance);
  CI95 = estimate + [-1 1] .* 1.96 .* SE;
else
  SE = NaN;
  CI95 = [NaN NaN];
end
result = struct('estimate', estimate, 'SE', SE, 'CI95', CI95);

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
