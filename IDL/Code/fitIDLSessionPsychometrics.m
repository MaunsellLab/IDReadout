function fit = fitIDLSessionPsychometrics(trialTable, varargin)
% fitIDLSessionPsychometrics  Jointly fit the enabled IDL allocation curves.
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

observed = sort(unique(round(double(trialTable.unsignedAllocPC), 6)));
assert(isequal(observed(:),[50;62.5;75;100]) || isequal(observed(:),[50;100]), ...
  'fitIDLSessionPsychometrics:UnexpectedAllocations', ...
  'Expected the four original allocations or 100:0/50:50; found %s.',mat2str(observed'));
allocationPC = sort(observed(:),'descend');
if ismember('splitUnsplitRatio',trialTable.Properties.VariableNames)
  two = numel(allocationPC)==2;
  ratio = double(trialTable.splitUnsplitRatio);
  assert(all(isfinite(ratio) & ratio>0) && ...
    all((abs(ratio-1)>1e-9)==two), ...
    'fitIDLSessionPsychometrics:Ratio','Allocation style and trial ratios disagree.');
end

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
    'fitIDLSessionPsychometrics:LapseBound', ...
    'The lapse upper bound must leave target performance below the asymptote.');
  initialLapse = min(max(0.05, lapseBounds(1) + eps), lapseBounds(2) - eps);
else
  fixedLapse = double(opts.FixedLapse);
  assert(targetPerformance < 1 - fixedLapse, ...
    'fitIDLSessionPsychometrics:FixedLapse', ...
    'Target performance must be below the fixed-lapse asymptote.');
  initialLapse = fixedLapse;
end

nAllocations = numel(allocationPC);
threshold0 = nan(nAllocations, 1);
for k = 1:nAllocations
  idx = abs(allocation - allocationPC(k)) < 1e-6;
  initialFit = fitIDQSessionWeibullFixedBeta( ...
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
  predictedP(idx) = idqWeibullP(double(trialTable.totalStepCohPC(idx)), ...
    alpha(k), betaWeibull, lapse);
end

fit = struct();
fit.model = 'Enabled allocation-specific thresholds; shared beta and lapse';
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
  p(idx) = idqWeibullP(coh(idx), alpha(k), betaWeibull, lapse);
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
alpha = idqWeibullAlphaForThreshold(targetPerformance, threshold, betaWeibull, lapse);
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
