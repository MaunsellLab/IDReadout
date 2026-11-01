function cfg = configureAdaptationExponent(cfg,exponent,referencePC)
%CONFIGUREADAPTATIONEXPONENT Select power-weighted prestep adaptation.
if nargin<3||isempty(referencePC)
  if isfield(cfg,'adaptationReferencePC'), referencePC=cfg.adaptationReferencePC;
  else, referencePC=40; end
end
exponent=double(exponent); referencePC=double(referencePC);
assert(isscalar(exponent)&&isfinite(exponent)&&exponent>=0, ...
  'AdaptationExponent must be a finite nonnegative scalar.');
assert(isscalar(referencePC)&&isfinite(referencePC)&&referencePC>0, ...
  'AdaptationReferencePC must be a positive finite scalar.');
cfg.adaptationExponent=exponent;
cfg.adaptationReferencePC=referencePC;
if ~isfield(cfg,'decisionRule'), cfg.decisionRule="absolute-step-max"; end
cfg=maxpoolfit.configureDecisionRule(cfg,cfg.decisionRule);
end
