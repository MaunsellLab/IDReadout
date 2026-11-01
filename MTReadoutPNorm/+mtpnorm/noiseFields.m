function Z = noiseFields(cfg,nDraws)
%NOISEFIELDS Four independent smooth directional Gaussian integration sets.
% Dimensions: direction x draw x RDK x epoch. Reusing these quadrature
% draws across trials/parameters stabilizes integration; it does not impose
% a shared physical noise realization on different trials.
assert(mod(nDraws,2)==0 && nDraws>=2,'Use an even number of draws.');
theta=deg2rad((0:cfg.directionSpacingDeg:360-cfg.directionSpacingDeg)');
ell=deg2rad(cfg.noiseCorrLengthDeg);
assert(ell>0,'Noise correlation length must be positive.');
K=exp(-2*sin((theta-theta')/2).^2/ell^2);
L=chol((K+K')/2+1e-10*eye(numel(theta)),'lower');
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(cfg.seed,'twister'); Z=zeros(numel(theta),nDraws,2,2);
for j=1:2
  for e=1:2
    z=L*randn(numel(theta),nDraws/2);
    Z(:,:,j,e)=[z -z];
  end
end
end
