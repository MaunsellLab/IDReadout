function D = decisionResponse(stepResponse,adaptationTemplate,preResponse,adaptation,cfg)
%DECISIONRESPONSE Construct directional evidence before hard-max pooling.
assert(isequal(size(stepResponse),size(adaptationTemplate),size(preResponse)), ...
  'Step, prestep, and adaptation response arrays must have matching sizes.');
if ~isfield(cfg,'decisionRule')
  rule="absolute-step-max";
else
  rule=string(cfg.decisionRule);
end
adaptationDrive=maxpoolfit.adaptationDrive(adaptationTemplate,cfg);
responseGain=maxpoolfit.adaptationGain(preResponse,cfg);
switch rule
  case "absolute-step-max"
    adaptedBaseline=preResponse-adaptation.*adaptationDrive;
    D=adaptedBaseline+responseGain.*(stepResponse-preResponse);
  case "directional-change-max"
    D=-adaptation.*adaptationDrive+ ...
      responseGain.*(stepResponse-preResponse);
  otherwise
    error('Unknown decision rule %s.',rule);
end
end
