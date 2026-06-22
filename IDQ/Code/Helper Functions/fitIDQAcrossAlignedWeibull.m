function fit = fitIDQAcrossAlignedWeibull(alignedCoh, correct, targetPerformance, lapseBounds)
% fitIDQAcrossAlignedWeibull
%
% Fit an across-session Weibull to threshold-aligned coherence values.
%
% Model:
%   P(c) = 0.5 + (0.5 - lapse) * (1 - exp(-(c / alpha)^beta))
%
% Free parameters:
%   threshold at targetPerformance, in aligned units
%   betaWeibull
%   lapse
%
% Inputs:
%   alignedCoh          threshold-aligned coherence values
%   correct             logical or 0/1 correctness
%   targetPerformance   e.g. 0.75
%   lapseBounds         e.g. [0 0.05]
%
% Output:
%   fit struct

alignedCoh = double(alignedCoh(:));
correct = double(correct(:));

valid = isfinite(alignedCoh) & isfinite(correct);
alignedCoh = alignedCoh(valid);
correct = correct(valid);

fit = struct();
fit.threshold = NaN;
fit.thresholdPerformance = targetPerformance;
fit.alpha = NaN;
fit.betaWeibull = NaN;
fit.lapse = NaN;
fit.nTrials = numel(correct);
fit.nCorrect = sum(correct == 1);
fit.negLogLikelihood = NaN;
fit.exitflag = NaN;
fit.hessian = [];
fit.paramNames = {'logThreshold', 'logBeta', 'lapseRaw'};

if numel(unique(alignedCoh)) < 2 || numel(unique(correct)) < 2
    return
end

positiveCoh = alignedCoh(alignedCoh > 0);
if isempty(positiveCoh)
    return
end

% Initial guesses. After alignment, threshold should be near 1.
threshold0 = 1;
beta0 = 3;
lapse0 = mean(lapseBounds);

theta0 = [
    log(threshold0)
    log(beta0)
    logitBounded(lapse0, lapseBounds)
    ];

objective = @(theta) negLogLikelihood(theta, alignedCoh, correct, ...
    targetPerformance, lapseBounds);

opts = optimset( ...
    'Display', 'off', ...
    'MaxFunEvals', 5000, ...
    'MaxIter', 5000);

[thetaHat, nll, exitflag] = fminsearch(objective, theta0, opts);

[threshold, betaWeibull, lapse, alpha] = unpackTheta( ...
    thetaHat, targetPerformance, lapseBounds);

fit.threshold = threshold;
fit.alpha = alpha;
fit.betaWeibull = betaWeibull;
fit.lapse = lapse;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.thetaHat = thetaHat;

% Optional finite-difference Hessian for later diagnostics.
fit.hessian = finiteDifferenceHessian(objective, thetaHat);

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihood(theta, coh, correct, targetPerformance, lapseBounds)

[~, betaWeibull, lapse, alpha] = unpackTheta(theta, targetPerformance, lapseBounds);

p = idqWeibullP(coh, alpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);

nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end

%% ------------------------------------------------------------------------
function [threshold, betaWeibull, lapse, alpha] = unpackTheta(theta, targetPerformance, lapseBounds)

threshold = exp(theta(1));
betaWeibull = exp(theta(2));
lapse = invLogitBounded(theta(3), lapseBounds);

alpha = idqWeibullAlphaForThreshold( ...
    targetPerformance, threshold, betaWeibull, lapse);

end

%% ------------------------------------------------------------------------
function y = logitBounded(x, bounds)

lo = bounds(1);
hi = bounds(2);

x = min(max(x, lo + eps), hi - eps);
u = (x - lo) ./ (hi - lo);

y = log(u ./ (1 - u));

end

%% ------------------------------------------------------------------------
function x = invLogitBounded(y, bounds)

lo = bounds(1);
hi = bounds(2);

u = 1 ./ (1 + exp(-y));
x = lo + (hi - lo) .* u;

end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(fun, theta)

theta = theta(:);
n = numel(theta);
H = nan(n, n);

f0 = fun(theta);
h = 1e-4;

for i = 1:n
    ei = zeros(n, 1);
    ei(i) = h;

    for j = i:n
        ej = zeros(n, 1);
        ej(j) = h;

        if i == j
            fp = fun(theta + ei);
            fm = fun(theta - ei);
            H(i, i) = (fp - 2 * f0 + fm) / h^2;
        else
            fpp = fun(theta + ei + ej);
            fpm = fun(theta + ei - ej);
            fmp = fun(theta - ei + ej);
            fmm = fun(theta - ei - ej);
            H(i, j) = (fpp - fpm - fmp + fmm) / (4 * h^2);
            H(j, i) = H(i, j);
        end
    end
end

end