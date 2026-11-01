function draws = noiseDraws(cfg, nDraws, seedOffset)
%NOISEDRAWS Standard correlated circular Gaussian-field integration nodes.
if nargin < 2 || isempty(nDraws), nDraws = cfg.nNoiseDraws; end
if nargin < 3, seedOffset = 0; end
if cfg.useAntitheticDraws
  assert(mod(nDraws,2)==0, 'Antithetic draw count must be even.');
  nBase = nDraws/2;
else
  nBase = nDraws;
end

theta = deg2rad(cfg.directionGridDeg(:));
ell = deg2rad(cfg.noiseCorrLengthDeg);
delta = theta-theta';
K = exp(-2 .* sin(delta./2).^2 ./ ell.^2);
K = (K+K')/2;
[L,p] = chol(K + 1e-10*eye(size(K)),'lower');
assert(p==0, 'Internal-noise covariance is not positive definite.');

old = rng;
cleanup = onCleanup(@()rng(old)); %#ok<NASGU>
rng(cfg.noiseSeed + seedOffset,'twister');
z = L*randn(numel(theta),nBase);
if cfg.useAntitheticDraws, z=[z -z]; end
draws = z;
end

