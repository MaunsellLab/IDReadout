function draws = independentNoiseDraws(cfg,nDraws,seedOffset)
%INDEPENDENTNOISEDRAWS Independent mechanism-noise integration nodes.
if nargin<2||isempty(nDraws), nDraws=cfg.nNoiseDraws; end
if nargin<3, seedOffset=0; end
if cfg.useAntitheticDraws
  assert(mod(nDraws,2)==0,'Antithetic draw count must be even.');
  nBase=nDraws/2;
else
  nBase=nDraws;
end
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(cfg.noiseSeed+seedOffset,'twister');
z=randn(numel(cfg.directionGridDeg),nBase);
if cfg.useAntitheticDraws, z=[z -z]; end
draws=z;
end
