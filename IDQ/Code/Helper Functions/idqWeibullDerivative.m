function dpdc = idqWeibullDerivative(coh, alpha, betaWeibull, lapse)
% idqWeibullDerivative
%
% Derivative of the IDQ Weibull psychometric function with respect to
% coherence.
%
%   p = 0.5 + (0.5 - lapse) * (1 - exp(-(coh / alpha)^beta))
%
%   dp/dc = (0.5 - lapse) * exp(-(coh/alpha)^beta) ...
%           * beta * coh^(beta - 1) / alpha^beta
%
% Inputs:
%   coh          coherence values, percent coherence units
%   alpha        Weibull scale parameter, same units as coh
%   betaWeibull  Weibull shape parameter
%   lapse        lapse rate
%
% Output:
%   dpdc         derivative in probability-correct per coherence unit

coh = double(coh);
alpha = double(alpha);
betaWeibull = double(betaWeibull);
lapse = double(lapse);

z = (coh ./ alpha) .^ betaWeibull;

dpdc = (0.5 - lapse) .* exp(-z) .* ...
  betaWeibull .* coh .^ (betaWeibull - 1) ./ alpha .^ betaWeibull;

% Handle the exact zero-coherence case explicitly.
if betaWeibull > 1
  dpdc(coh == 0) = 0;
elseif betaWeibull == 1
  dpdc(coh == 0) = (0.5 - lapse) .* betaWeibull ./ alpha;
end

end