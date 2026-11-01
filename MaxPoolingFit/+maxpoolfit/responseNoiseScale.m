function [scale,fullResponse,referenceResponse] = responseNoiseScale(B,cfg)
%RESPONSENOISESCALE Scale independent variance with full sensory response.
% B is the zero-centered fixed-density step response. The response to the
% isotropic random-dot pool supplies the physically defined positive offset.
p=0;
if isfield(cfg,'noiseVarianceExponent')
  p=double(cfg.noiseVarianceExponent);
end
assert(isscalar(p)&&isfinite(p)&&p>=0, ...
  'cfg.noiseVarianceExponent must be a finite nonnegative scalar.');
referenceResponse=100*maxpoolfit.randomDirectionMean(cfg.sigmaMTDeg);
fullResponse=referenceResponse+B;
assert(all(isfinite(fullResponse),'all'), ...
  'Full sensory responses must be finite.');
if p==0
  % Exact homoscedastic nesting: the response proxy is diagnostically
  % irrelevant and must not impose additional stimulus restrictions.
  scale=ones(size(B),'like',B);
else
  % The coherence-component representation is linear and can extrapolate
  % below zero (for example, nominal components can sum above 100%). A
  % firing rate cannot. Rectification is therefore the appropriate
  % Poisson-like response proxy; eps only avoids an exact numerical zero.
  fullResponse=max(fullResponse,eps(referenceResponse));
  % p is the variance-response exponent, so SD scales with p/2.
  scale=(fullResponse./referenceResponse).^(p/2);
end
end
