function fit = fitIDQSessionWeibullFixedBeta(stepCoh, correct, betaWeibull, lapse, targetPerformance)
% fitIDQSessionWeibullFixedBeta
%
% Fit a single IDQ session Weibull with fixed beta and lapse.
% The only fitted parameter is threshold at targetPerformance.

stepCoh = double(stepCoh(:));
correct = double(correct(:));

valid = isfinite(stepCoh) & isfinite(correct);
stepCoh = stepCoh(valid);
correct = correct(valid);

fit = struct();
fit.threshold = NaN;
fit.thresholdPerformance = targetPerformance;
fit.alpha = NaN;
fit.betaWeibull = betaWeibull;
fit.lapse = lapse;
fit.nTrials = numel(correct);
fit.nCorrect = sum(correct == 1);
fit.negLogLikelihood = NaN;
fit.exitflag = NaN;

if numel(unique(stepCoh)) < 2 || numel(unique(correct)) < 2
  return
end

positiveCoh = stepCoh(stepCoh > 0);
if isempty(positiveCoh)
  return
end

lowerThreshold = max(min(positiveCoh) / 10, eps);
upperThreshold = max(positiveCoh) * 10;

lb = log(lowerThreshold);
ub = log(upperThreshold);

objective = @(theta) negLogLikelihood(theta, stepCoh, correct, ...
  betaWeibull, lapse, targetPerformance);

opts = optimset( ...
  'Display', 'off', ...
  'MaxFunEvals', 2000, ...
  'MaxIter', 2000);

[thetaHat, nll, exitflag] = fminbnd(objective, lb, ub, opts);

threshold = exp(thetaHat);
alpha = idqWeibullAlphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);

fit.threshold = threshold;
fit.alpha = alpha;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihood(theta, stepCoh, correct, betaWeibull, lapse, targetPerformance)

threshold = exp(theta);
alpha = idqWeibullAlphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);

p = idqWeibullP(stepCoh, alpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);

nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end