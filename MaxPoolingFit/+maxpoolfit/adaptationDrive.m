function drive = adaptationDrive(prestepResponse,cfg)
%ADAPTATIONDRIVE Apply nested power weighting to directional prestep drive.
%
% drive = Ppre .* (Ppre/P0).^q.  q=0 exactly recovers the linear
% prestep-adaptation model. P0 fixes the interpretation of the shared
% adaptation coefficient but does not alter the direction-specific form.

q=0;
if isfield(cfg,'adaptationExponent'), q=double(cfg.adaptationExponent); end
P0=40;
if isfield(cfg,'adaptationReferencePC'), P0=double(cfg.adaptationReferencePC); end
assert(isscalar(q)&&isfinite(q)&&q>=0, ...
  'cfg.adaptationExponent must be a finite nonnegative scalar.');
assert(isscalar(P0)&&isfinite(P0)&&P0>0, ...
  'cfg.adaptationReferencePC must be a positive finite scalar.');
assert(all(isfinite(prestepResponse),'all'), ...
  'Directional prestep responses must be finite.');
if q==0
  drive=prestepResponse;
else
  assert(all(prestepResponse>=0,'all'), ...
    ['Power-weighted adaptation requires a nonnegative full response. ' ...
     'It is not defined for zero-centered redistribution responses.']);
  drive=prestepResponse.*(prestepResponse./P0).^q;
end
end
