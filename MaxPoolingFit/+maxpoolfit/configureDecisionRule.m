function cfg = configureDecisionRule(cfg,rule)
%CONFIGUREDECISIONRULE Select the between-RDK evidence computation.
rule=string(rule);
assert(isscalar(rule) && ismember(rule, ...
  ["absolute-step-max" "directional-change-max"]), ...
  'DecisionRule must be "absolute-step-max" or "directional-change-max".');
cfg.decisionRule=rule;
redistribution=isfield(cfg,'coherenceResponseMode') && ...
  string(cfg.coherenceResponseMode)=="random-dot-redistribution";
powerModel=isfield(cfg,'adaptationExponent') && ...
  abs(double(cfg.adaptationExponent))>1e-12;
gainModel=isfield(cfg,'adaptationGainStrength') && ...
  abs(double(cfg.adaptationGainStrength))>1e-12;
if rule=="absolute-step-max" && ~powerModel && ~gainModel
  cfg.modelVersion="DenseHardMax-v5-shared-prestep-adaptation";
elseif rule=="absolute-step-max" && powerModel && ~gainModel
  cfg.modelVersion="DenseHardMax-v8-power-prestep-adaptation";
elseif rule=="absolute-step-max" && ~powerModel
  cfg.modelVersion="DenseHardMax-v14-increment-only-adaptation-gain";
elseif rule=="absolute-step-max"
  cfg.modelVersion="DenseHardMax-v16-power-and-increment-gain";
elseif ~powerModel && ~gainModel
  cfg.modelVersion="DenseHardMax-v6-directional-change";
elseif powerModel && ~gainModel
  cfg.modelVersion="DenseHardMax-v9-directional-change-power-adaptation";
elseif ~powerModel
  cfg.modelVersion="DenseHardMax-v15-directional-increment-adaptation-gain";
else
  cfg.modelVersion="DenseHardMax-v17-directional-power-and-increment-gain";
end
if redistribution
  cfg.modelVersion=cfg.modelVersion+"-random-dot-redistribution";
end
end
