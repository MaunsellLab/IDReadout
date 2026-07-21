function signedFits = fitIDQSignedNoiseAsymmetry(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDQSignedNoiseAsymmetry
%
% Test whether change-side drift and pooled non-drift noise have different
% slopes for positive and negative predictor values. Predictor signs remain
% physical coherence signs.
%
% The null model is the existing reduced main-effects model:
%   cEff = cStep + gD*nD + gN*(nPlus+nMinus)
%
% Signed asymmetry is represented by centered absolute-value predictors:
%   aN = |nPlus| + |nMinus|, centered within session
%   aD = |nD|, centered within session
%
% With hN*aN in the model, the non-drift slopes are:
%   negative noise: gN - hN
%   positive noise: gN + hN
% Centering aN and aD within session prevents asymmetry terms from changing
% a session's mean effective coherence.

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "noisePredDir1", "noisePredDir2", "noisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQSignedNoiseAsymmetry:MissingVariables', ...
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
  error('fitIDQSignedNoiseAsymmetry:MissingSessionThreshold', ...
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

aN = centerWithinSession(abs(nPlus) + abs(nMinus), sessionIndex);
aD = centerWithinSession(abs(nD), sessionIndex);

common = struct();
common.stepCoh = stepCoh;
common.correct = correct;
common.sessionAlpha = sessionAlpha;
common.betaWeibull = betaWeibull;
common.lapse = lapse;

Xmain = [nD, nN];
main = fitGainModel('main', ["g_D", "g_N"], Xmain, [1; 0], common);
nonDriftAsymmetry = fitGainModel('nonDriftAsymmetry', ...
  ["g_D", "g_N", "h_N"], [Xmain, aN], [main.gain; 0], common);
driftAsymmetry = fitGainModel('driftAsymmetry', ...
  ["g_D", "g_N", "h_D"], [Xmain, aD], [main.gain; 0], common);
bothAsymmetries = fitGainModel('bothAsymmetries', ...
  ["g_D", "g_N", "h_N", "h_D"], [Xmain, aN, aD], ...
  [main.gain; 0; 0], common);

nonDriftAsymmetry = addNestedComparison(nonDriftAsymmetry, main);
driftAsymmetry = addNestedComparison(driftAsymmetry, main);
bothAsymmetries = addNestedComparison(bothAsymmetries, main);

signedFits = struct();
signedFits.main = main;
signedFits.nonDriftAsymmetry = nonDriftAsymmetry;
signedFits.driftAsymmetry = driftAsymmetry;
signedFits.bothAsymmetries = bothAsymmetries;
signedFits.nonDriftSignedSlopes = makeSignedSlopes(nonDriftAsymmetry, 2, 3, 'Non-Drift');
signedFits.driftSignedSlopes = makeSignedSlopes(driftAsymmetry, 1, 3, 'Drift');
signedFits.jointNonDriftSignedSlopes = makeSignedSlopes(bothAsymmetries, 2, 3, 'Non-Drift');
signedFits.jointDriftSignedSlopes = makeSignedSlopes(bothAsymmetries, 1, 4, 'Drift');

signedFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nD, nPlus, nMinus, nN, aN, aD, main.pPred, ...
  nonDriftAsymmetry.pPred, driftAsymmetry.pPred, bothAsymmetries.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nD', 'nPlus', 'nMinus', 'nN', 'aN', 'aD', 'pMain', ...
  'pNonDriftAsymmetry', 'pDriftAsymmetry', 'pBothAsymmetries'});
signedFits.nTrials = nTrials;
signedFits.betaWeibull = betaWeibull;
signedFits.lapse = lapse;
signedFits.targetPerformance = targetPerformance;
signedFits.predictorDefinition = ...
  ['Absolute-value asymmetry predictors are centered within session. ' ...
   'For each stream, gPositive = g + h and gNegative = g - h.'];

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
function slopes = makeSignedSlopes(fit, gainIndex, asymmetryIndex, streamName)

Lnegative = zeros(1, fit.nParameters);
Lpositive = zeros(1, fit.nParameters);
Lnegative(gainIndex) = 1;
Lnegative(asymmetryIndex) = -1;
Lpositive(gainIndex) = 1;
Lpositive(asymmetryIndex) = 1;

negative = makeLinearCombination(fit, Lnegative);
positive = makeLinearCombination(fit, Lpositive);

signName = ["Negative"; "Positive"];
gain = [negative.estimate; positive.estimate];
SE = [negative.SE; positive.SE];
CI95Low = [negative.CI95(1); positive.CI95(1)];
CI95High = [negative.CI95(2); positive.CI95(2)];
slopes = table(signName, gain, SE, CI95Low, CI95High);
slopes.Properties.Description = sprintf('%s signed-noise slopes', streamName);

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
