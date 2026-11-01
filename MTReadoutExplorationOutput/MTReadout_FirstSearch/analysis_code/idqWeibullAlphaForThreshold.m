function alpha = idqWeibullAlphaForThreshold(pTarget, threshold, betaWeibull, lapse)
% idqWeibullAlphaForThreshold
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
  error('idqWeibullAlphaForThreshold:InvalidTarget', ...
    'pTarget must lie between 0.5 and 1 - lapse.');
end

alpha = threshold ./ (-log(1 - q)) .^ (1 ./ betaWeibull);

end