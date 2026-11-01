function fit = fitIDJSessionPsychometrics(trialTable, varargin)
% fitIDJSessionPsychometrics  Jointly fit the four IDJ allocation curves.
%
% Fits allocation-specific thresholds with one shared Weibull beta and one
% shared lapse. Pass 'FixedLapse', 0 for a zero-lapse fit or leave it empty
% to estimate a common lapse within LapseBounds.

p = inputParser;
p.addParameter('TargetPerformance', 0.75, ...
  @(x) isnumeric(x) && isscalar(x) && x > 0.5 && x < 1);
p.addParameter('FixedLapse', [], ...
  @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
p.addParameter('LapseBounds', [0 0.15], ...
  @(x) isnumeric(x) && numel(x) == 2 && x(1) >= 0 && x(2) > x(1));
p.parse(varargin{:});
opts = p.Results;

assert(ismember('topUpTrial',trialTable.Properties.VariableNames) && ...
  all(double(trialTable.topUpTrial)==0), ...
  'fitIDJSessionPsychometrics:TopUpTrials', ...
  'Input trials must have topUpTrial=false. Rebuild the compact session analysis.');

expected = [100; 70; 50; 100/3];
observed = unique(double(trialTable.unsignedAllocPC));
assert(all(any(abs(observed(:) - expected') < 1e-6, 2)), ...
  'fitIDJSessionPsychometrics:UnexpectedAllocations', 'Unexpected IDJ allocation.');
allocationPC = expected(any(abs(expected - observed(:)') < 1e-6, 2));
assert(~isempty(allocationPC), 'fitIDJSessionPsychometrics:NoTrials', ...
  'No allocation data are available.');

coh = double(trialTable.totalStepCohPC(:));
correct = double(trialTable.correct(:));
allocation = double(trialTable.unsignedAllocPC(:));
valid = isfinite(coh) & isfinite(correct) & isfinite(allocation);
coh = coh(valid);
correct = correct(valid);
allocation = allocation(valid);

targetPerformance = opts.TargetPerformance;
lapseBounds = double(opts.LapseBounds(:)');
if isempty(opts.FixedLapse)
  assert(targetPerformance < 1 - lapseBounds(2), ...
    'fitIDJSessionPsychometrics:LapseBound', ...
    'The lapse upper bound must leave target performance below the asymptote.');
  initialLapse = min(max(0.05, lapseBounds(1) + eps), lapseBounds(2) - eps);
else
  fixedLapse = double(opts.FixedLapse);
  assert(targetPerformance < 1 - fixedLapse, ...
    'fitIDJSessionPsychometrics:FixedLapse', ...
    'Target performance must be below the fixed-lapse asymptote.');
  initialLapse = fixedLapse;
end

nAllocations = numel(allocationPC);
threshold0 = nan(nAllocations, 1);
for k = 1:nAllocations
  idx = abs(allocation - allocationPC(k)) < 1e-6;
  initialFit = initialFixedBetaFit( ...
    coh(idx), correct(idx), 2, initialLapse, targetPerformance);
  threshold0(k) = initialFit.threshold;
  if ~isfinite(threshold0(k))
    positiveCoh = coh(idx & coh > 0);
    threshold0(k) = median(positiveCoh);
  end
end

theta0 = [log(threshold0); log(2)];
if isempty(opts.FixedLapse)
  theta0(end + 1, 1) = logitBounded(initialLapse, lapseBounds);
end

objective = @(theta) negLogLikelihood(theta, coh, correct, allocation, ...
  allocationPC, targetPerformance, opts.FixedLapse, lapseBounds);
optimOptions = optimset('Display', 'off', 'MaxFunEvals', 10000, 'MaxIter', 10000);
[thetaHat, nll, exitflag] = fminsearch(objective, theta0, optimOptions);
[threshold, betaWeibull, lapse, alpha] = unpackTheta(thetaHat, nAllocations, ...
  targetPerformance, opts.FixedLapse, lapseBounds);

hessian = finiteDifferenceHessian(objective, thetaHat);
covariance = nan(numel(thetaHat));
if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;
  if all(eig(hessian) > 0)
    covariance = inv(hessian);
  end
end

nTrials = nan(nAllocations, 1);
nCorrect = nan(nAllocations, 1);
meanCorrect = nan(nAllocations, 1);
thresholdSE = nan(nAllocations, 1);
thresholdCILow = nan(nAllocations, 1);
thresholdCIHigh = nan(nAllocations, 1);
for k = 1:nAllocations
  idx = abs(allocation - allocationPC(k)) < 1e-6;
  nTrials(k) = sum(idx);
  nCorrect(k) = sum(correct(idx));
  meanCorrect(k) = mean(correct(idx));
  if all(isfinite(covariance), 'all') && covariance(k,k) > 0
    logSE = sqrt(covariance(k,k));
    thresholdSE(k) = threshold(k) * logSE;
    thresholdCILow(k) = threshold(k) * exp(-1.96 * logSE);
    thresholdCIHigh(k) = threshold(k) * exp(1.96 * logSE);
  end
end

thresholdTable = table(allocationPC, nTrials, nCorrect, meanCorrect, threshold, ...
  thresholdSE, thresholdCILow, thresholdCIHigh, alpha);

betaSE = NaN;
betaCI95 = [NaN NaN];
betaIndex = nAllocations + 1;
if all(isfinite(covariance), 'all') && covariance(betaIndex,betaIndex) > 0
  logSE = sqrt(covariance(betaIndex,betaIndex));
  betaSE = betaWeibull * logSE;
  betaCI95 = betaWeibull * exp([-1 1] * 1.96 * logSE);
end

lapseSE = 0;
lapseCI95 = [lapse lapse];
if isempty(opts.FixedLapse)
  lapseIndex = nAllocations + 2;
  lapseSE = NaN;
  lapseCI95 = [NaN NaN];
  if all(isfinite(covariance), 'all') && covariance(lapseIndex,lapseIndex) > 0
    u = (lapse - lapseBounds(1)) / diff(lapseBounds);
    derivative = diff(lapseBounds) * u * (1 - u);
    lapseSE = abs(derivative) * sqrt(covariance(lapseIndex,lapseIndex));
    lapseCI95 = min(max(lapse + [-1 1] * 1.96 * lapseSE, lapseBounds(1)), lapseBounds(2));
  end
end

predictedP = nan(height(trialTable), 1);
for k = 1:nAllocations
  idx = abs(double(trialTable.unsignedAllocPC) - allocationPC(k)) < 1e-6;
  predictedP(idx) = idjWeibullP(double(trialTable.totalStepCohPC(idx)), ...
    alpha(k), betaWeibull, lapse);
end

fit = struct();
fit.model = sprintf('%d allocation-specific thresholds; shared beta and lapse', nAllocations);
fit.targetPerformance = targetPerformance;
fit.fixedLapse = opts.FixedLapse;
fit.lapseBounds = lapseBounds;
fit.thresholdTable = thresholdTable;
fit.betaWeibull = betaWeibull;
fit.betaSE = betaSE;
fit.betaCI95 = betaCI95;
fit.lapse = lapse;
fit.lapseSE = lapseSE;
fit.lapseCI95 = lapseCI95;
fit.predictedP = predictedP;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.thetaHat = thetaHat;
fit.hessian = hessian;
fit.covariance = covariance;
fit.nTrials = numel(correct);
fit.nCorrect = sum(correct);
end

%% ------------------------------------------------------------------------
function nll = negLogLikelihood(theta, coh, correct, allocation, allocationPC, ...
    targetPerformance, fixedLapse, lapseBounds)

[~, betaWeibull, lapse, alpha] = unpackTheta(theta, numel(allocationPC), ...
  targetPerformance, fixedLapse, lapseBounds);
p = nan(size(correct));
for k = 1:numel(allocationPC)
  idx = abs(allocation - allocationPC(k)) < 1e-6;
  p(idx) = idjWeibullP(coh(idx), alpha(k), betaWeibull, lapse);
end
p = min(max(p, eps), 1 - eps);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));
end

%% ------------------------------------------------------------------------
function [threshold, betaWeibull, lapse, alpha] = unpackTheta(theta, ...
    nAllocations, targetPerformance, fixedLapse, lapseBounds)

threshold = exp(theta(1:nAllocations));
betaWeibull = exp(theta(nAllocations + 1));
if isempty(fixedLapse)
  lapse = invLogitBounded(theta(nAllocations + 2), lapseBounds);
else
  lapse = double(fixedLapse);
end
alpha = alphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);
end

%% ------------------------------------------------------------------------
function transformed = logitBounded(value, bounds)

u = (value - bounds(1)) / diff(bounds);
transformed = log(u / (1 - u));
end

%% ------------------------------------------------------------------------
function value = invLogitBounded(transformed, bounds)

u = 1 / (1 + exp(-transformed));
value = bounds(1) + diff(bounds) * u;
end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(fun, theta)

theta = theta(:);
n = numel(theta);
H = nan(n);
f0 = fun(theta);
h = 1e-4 .* max(1, abs(theta));
for i = 1:n
  ei = zeros(n,1);
  ei(i) = h(i);
  H(i,i) = (fun(theta + ei) - 2 * f0 + fun(theta - ei)) / h(i)^2;
  for j = i + 1:n
    ej = zeros(n,1);
    ej(j) = h(j);
    H(i,j) = (fun(theta + ei + ej) - fun(theta + ei - ej) - ...
      fun(theta - ei + ej) + fun(theta - ei - ej)) / (4 * h(i) * h(j));
    H(j,i) = H(i,j);
  end
end
end

function fit = initialFixedBetaFit(stepCoh, correct, betaWeibull, lapse, targetPerformance)
% initialFixedBetaFit
%
% Fit a single IDJ session Weibull with fixed beta and lapse.
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

objective = @(theta) initialNegLogLikelihood(theta, stepCoh, correct, ...
  betaWeibull, lapse, targetPerformance);

opts = optimset( ...
  'Display', 'off', ...
  'MaxFunEvals', 2000, ...
  'MaxIter', 2000);

[thetaHat, nll, exitflag] = fminbnd(objective, lb, ub, opts);

threshold = exp(thetaHat);
alpha = alphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);

fit.threshold = threshold;
fit.alpha = alpha;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;

end

%% ------------------------------------------------------------------------
function nll = initialNegLogLikelihood(theta, stepCoh, correct, betaWeibull, lapse, targetPerformance)

threshold = exp(theta);
alpha = alphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);

p = idjWeibullP(stepCoh, alpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);

nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end
function alpha = alphaForThreshold(pTarget, threshold, betaWeibull, lapse)
% alphaForThreshold
%
% Return the Weibull alpha/scale parameter implied by a threshold defined as
% the coherence giving pTarget.
%
%   p = 0.5 + (0.5 - lapse) * (1 - exp(-(coh / alpha)^beta))
%
% Inputs:
%   pTarget      target probability correct, e.g. 0.75
%   threshold    coherence giving pTarget
%   betaWeibull  Weibull shape parameter
%   lapse        lapse rate
%
% Output:
%   alpha        Weibull scale parameter
%
% If:
%
%   threshold = alpha * (-log(1 - q))^(1 / beta)
%
% then:
%
%   alpha = threshold / (-log(1 - q))^(1 / beta)

pTarget = double(pTarget);
threshold = double(threshold);
betaWeibull = double(betaWeibull);
lapse = double(lapse);

q = (pTarget - 0.5) ./ (0.5 - lapse);

if any(q <= 0 | q >= 1, 'all')
  error('alphaForThreshold:InvalidTarget', ...
    'pTarget must lie between 0.5 and 1 - lapse.');
end

alpha = threshold ./ (-log(1 - q)) .^ (1 ./ betaWeibull);

end
