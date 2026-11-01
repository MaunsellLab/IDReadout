function gainFits = fitIDSNoiseGain(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDSNoiseGain  Fit IDS drift-relative coherence-noise gains.
%
% Predictor order is invariant across trials:
%   Drift, Drift+Probe (Probe0), Drift-Probe (Probe1).

idx = trialTable.hasStepNoise;
sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
X = double([trialTable.noisePredDrift(idx), trialTable.noisePredPlus(idx), ...
  trialTable.noisePredMinus(idx)]);
valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & all(isfinite(X), 2);
sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
X = X(valid, :);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;
sessionThreshold = nan(size(sessionIndex));
for k = 1:height(sessionFits)
  sessionThreshold(sessionIndex == sessionFits.sessionIndex(k)) = sessionFits.threshold(k);
end
if any(~isfinite(sessionThreshold))
  error('fitIDSNoiseGain:MissingSessionThreshold', ...
    'Could not assign a threshold to every noisy trial.');
end
sessionAlpha = idqWeibullAlphaForThreshold(targetPerformance, sessionThreshold, betaWeibull, lapse);

common = struct('stepCoh', stepCoh, 'correct', correct, ...
  'sessionAlpha', sessionAlpha, 'betaWeibull', betaWeibull, 'lapse', lapse);
gainFits = struct();
gainFits.combined = fitGainModel('combined', "All Streams", sum(X, 2), 1/3, common);
gainFits.driftRelative = fitGainModel('driftRelative', ...
  ["Drift", "+Probe", "-Probe"], X, [1; 0; 0], common);
gainFits.driftProbe = fitGainModel('driftProbe', ...
  ["Drift", "Probe"], [X(:,1), X(:,2) + X(:,3)], [1; 0], common);
gainFits.normalized = makeNormalizedGains(gainFits.driftRelative, gainFits.driftProbe);
end

function fit = fitGainModel(modelName, predictorNames, X, gain0, common)
nParameters = size(X, 2);
lb = -5 * ones(nParameters, 1);
ub = 5 * ones(nParameters, 1);
objective = @(gain) gainNLL(gain, X, common);
opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
  'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10, ...
  'MaxFunctionEvaluations', 5000, 'MaxIterations', 2000);
[gainHat, nll, exitflag] = fmincon(objective, gain0(:), [], [], [], [], lb, ub, [], opts);
gainHat = gainHat(:);
hessian = finiteDifferenceHessian(objective, gainHat);
covariance = nan(nParameters);
SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);
if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;
  if all(eig(hessian) > 0)
    covariance = inv(hessian);
    variance = diag(covariance);
    ok = isfinite(variance) & variance > 0;
    SE(ok) = sqrt(variance(ok));
    CI95(ok,:) = gainHat(ok) + [-1 1] .* (1.96 * SE(ok));
  end
end
fit = struct('model', modelName, 'predictorNames', predictorNames(:)', ...
  'gain', gainHat, 'SE', SE, 'CI95', CI95, 'covariance', covariance, ...
  'hessian', hessian, 'negLogLikelihood', nll, 'exitflag', exitflag, ...
  'nParameters', nParameters, 'nTrials', numel(common.correct), ...
  'nEffectiveCohClipped', sum(common.stepCoh + X * gainHat < 0));
end

function nll = gainNLL(gain, X, common)
effectiveCoh = max(common.stepCoh + X * gain(:), 0);
p = idqWeibullP(effectiveCoh, common.sessionAlpha, common.betaWeibull, common.lapse);
p = min(max(p, eps), 1 - eps);
nll = -sum(common.correct .* log(p) + (1-common.correct) .* log(1-p));
end

function N = makeNormalizedGains(relativeFit, pooledFit)
N = struct();
[N.pooledProbeOverDrift, N.pooledSE, N.pooledCI95] = ratioEstimate( ...
  pooledFit.gain, pooledFit.covariance, 2, 1);
[N.plusProbeOverDrift, N.plusSE, N.plusCI95] = ratioEstimate( ...
  relativeFit.gain, relativeFit.covariance, 2, 1);
[N.minusProbeOverDrift, N.minusSE, N.minusCI95] = ratioEstimate( ...
  relativeFit.gain, relativeFit.covariance, 3, 1);
N.method = 'Delta method using the joint observed-Hessian covariance';
end

function [ratio, se, ci] = ratioEstimate(gain, covariance, numerator, denominator)
ratio = gain(numerator) / gain(denominator);
se = NaN;
ci = [NaN NaN];
if gain(denominator) == 0 || any(~isfinite(covariance), 'all')
  return
end
gradient = zeros(numel(gain), 1);
gradient(numerator) = 1 / gain(denominator);
gradient(denominator) = -gain(numerator) / gain(denominator)^2;
variance = gradient' * covariance * gradient;
if isfinite(variance) && variance > 0
  se = sqrt(variance);
  ci = ratio + [-1 1] * 1.96 * se;
end
end

function H = finiteDifferenceHessian(fun, x)
x = x(:); n = numel(x); H = nan(n); h = 3e-4 .* max(1, abs(x)); f0 = fun(x);
for i = 1:n
  ei = zeros(n,1); ei(i) = h(i);
  H(i,i) = (fun(x+ei) - 2*f0 + fun(x-ei)) / h(i)^2;
  for j = i+1:n
    ej = zeros(n,1); ej(j) = h(j);
    H(i,j) = (fun(x+ei+ej)-fun(x+ei-ej)-fun(x-ei+ej)+fun(x-ei-ej)) / (4*h(i)*h(j));
    H(j,i) = H(i,j);
  end
end
end
