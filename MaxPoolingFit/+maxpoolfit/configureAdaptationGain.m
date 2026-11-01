function cfg = configureAdaptationGain(cfg,strength,referencePC)
%CONFIGUREADAPTATIONGAIN Select prestep-dependent exponential response gain.
if nargin<3||isempty(referencePC)
  if isfield(cfg,'adaptationReferencePC'), referencePC=cfg.adaptationReferencePC;
  else, referencePC=40; end
end
strength=double(strength); referencePC=double(referencePC);
assert(isscalar(strength)&&isfinite(strength)&&strength>=0, ...
  'AdaptationGainStrength must be a finite nonnegative scalar.');
assert(isscalar(referencePC)&&isfinite(referencePC)&&referencePC>0, ...
  'AdaptationReferencePC must be a positive finite scalar.');
cfg.adaptationGainStrength=strength;
cfg.adaptationReferencePC=referencePC;
if ~isfield(cfg,'decisionRule'), cfg.decisionRule="absolute-step-max"; end
cfg=maxpoolfit.configureDecisionRule(cfg,cfg.decisionRule);
end
