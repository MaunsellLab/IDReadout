function gain = adaptationGain(prestepResponse,cfg)
%ADAPTATIONGAIN Prestep-dependent multiplicative response gain.
%
% gain = exp(-b*Ppre/P0). b=0 exactly recovers unit gain.

b=0;
if isfield(cfg,'adaptationGainStrength')
  b=double(cfg.adaptationGainStrength);
end
P0=40;
if isfield(cfg,'adaptationReferencePC')
  P0=double(cfg.adaptationReferencePC);
end
assert(isscalar(b)&&isfinite(b)&&b>=0, ...
  'cfg.adaptationGainStrength must be a finite nonnegative scalar.');
assert(isscalar(P0)&&isfinite(P0)&&P0>0, ...
  'cfg.adaptationReferencePC must be a positive finite scalar.');
assert(all(isfinite(prestepResponse),'all'), ...
  'Directional prestep responses must be finite.');
if b==0
  gain=ones(size(prestepResponse),'like',prestepResponse);
else
  assert(all(prestepResponse>=0,'all'), ...
    ['Multiplicative adaptation gain requires a nonnegative full response. ' ...
     'It is not defined for zero-centered redistribution responses.']);
  gain=exp(-b.*prestepResponse./P0);
end
end
