function p = idjWeibullP(coh, alpha, betaWeibull, lapse)
% idjWeibullP  IDJ 2AFC Weibull; chance 0.5, upper asymptote 1-lapse.
z = (double(coh) ./ double(alpha)) .^ double(betaWeibull);
p = 0.5 + (0.5 - double(lapse)) .* (1 - exp(-z));
end
