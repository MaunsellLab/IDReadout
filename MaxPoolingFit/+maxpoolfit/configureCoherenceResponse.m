function cfg = configureCoherenceResponse(cfg,mode)
%CONFIGURECOHERENCERESPONSE Select the fixed-density RDK response contract.
mode=string(mode);
assert(isscalar(mode)&&ismember(mode, ...
  ["additive" "random-dot-redistribution"]), ...
  'CoherenceResponseMode must be "additive" or "random-dot-redistribution".');
cfg.coherenceResponseMode=mode;
if ~isfield(cfg,'decisionRule'), cfg.decisionRule="absolute-step-max"; end
cfg=maxpoolfit.configureDecisionRule(cfg,cfg.decisionRule);
end
