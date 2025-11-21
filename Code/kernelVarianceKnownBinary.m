function varK = kernelVarianceKnownBinary(X, choices)
% kernelVarianceKnownBinary
% Analytic per-bin variance of a mean-based kernel for known binary noise.
%
%   X:        T x N matrix of noise values (±a)
%   choices:  1 x N vector of trial labels (+1 / -1)
%
% Returns:
%   varK:     variance of all individual kernel bins

    % scalar amplitude of the binary noise
    a = mean(abs(X(:,1)));  % assumes all trials use the same amplitude

    % Count trials in each choice group (same for all time bins if no censoring)
    n_plus  = sum((choices ==  1));
    n_minus = sum((choices == -1));

    % Known stimulus variance for ±a, p = 0.5
    sigma2 = a^2;

    % Variance of kernel at each time bin: sigma^2 * (1/n_plus + 1/n_minus)
    varK = sigma2 * (1/n_plus + 1/n_minus);
end
