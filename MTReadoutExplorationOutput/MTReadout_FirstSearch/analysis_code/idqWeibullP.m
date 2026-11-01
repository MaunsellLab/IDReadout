function p = idqWeibullP(coh, alpha, betaWeibull, lapse)
% idqWeibullP
%
% IDQ 2AFC Weibull psychometric function.
%
%   p = 0.5 + (0.5 - lapse) * (1 - exp(-(coh / alpha)^beta))
%
% Inputs:
%   coh          coherence values, percent coherence units
%   alpha        Weibull scale parameter, same units as coh
%   betaWeibull  Weibull shape parameter
%   lapse        lapse rate; upper asymptote is 1 - lapse
%
% Output:
%   p            predicted probability correct
%
% Notes:
%   The function is anchored at p = 0.5 when coh = 0.
%   For betaWeibull > 1, slope at coh = 0 is zero.

coh = double(coh);
alpha = double(alpha);
betaWeibull = double(betaWeibull);
lapse = double(lapse);

z = (coh ./ alpha) .^ betaWeibull;

p = 0.5 + (0.5 - lapse) .* (1 - exp(-z));

end