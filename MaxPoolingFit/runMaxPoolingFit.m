function result = runMaxPoolingFit(cachePath, outputFolder, varargin)
%RUNMAXPOOLINGFIT Development and final domain-noise fits and predictions.
%
% result = runMaxPoolingFit(cachePath, outputFolder)

p=inputParser;
defaultCfg=maxpoolfit.defaultConfig();
p.addParameter('DevelopmentDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('FinalDraws',256,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('Adaptation',defaultCfg.adaptation, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
p.addParameter('AdaptationExponent',defaultCfg.adaptationExponent, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.addParameter('AdaptationGainStrength',defaultCfg.adaptationGainStrength, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.addParameter('AdaptationReferencePC',defaultCfg.adaptationReferencePC, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.addParameter('AdaptationProfile',[],@(x)isempty(x)||isstruct(x));
p.addParameter('DecisionRule',defaultCfg.decisionRule, ...
  @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end

cache=maxpoolfit.loadCache(cachePath);
cfg=defaultCfg;
cfg=maxpoolfit.configureDecisionRule(cfg,p.Results.DecisionRule);
cfg=maxpoolfit.configureAdaptationExponent(cfg,p.Results.AdaptationExponent, ...
  p.Results.AdaptationReferencePC);
cfg=maxpoolfit.configureAdaptationGain(cfg,p.Results.AdaptationGainStrength, ...
  p.Results.AdaptationReferencePC);
cfg.adaptation=double(p.Results.Adaptation);
if abs(cfg.adaptation-defaultCfg.adaptation)>1e-12
  cfg.adaptationSelection="explicit runMaxPoolingFit Adaptation input";
end
cfg.nNoiseDraws=p.Results.DevelopmentDraws;
developmentFit=maxpoolfit.fitModel(cache,cfg,'NDraws',p.Results.DevelopmentDraws);

cfg.nNoiseDraws=p.Results.FinalDraws;
cfg.finalNoiseDraws=p.Results.FinalDraws;
finalFit=maxpoolfit.fitModel(cache,cfg,'NDraws',p.Results.FinalDraws);
prediction=maxpoolfit.predict(cache,finalFit,cfg, ...
  'TrialMask',finalFit.trialMask,'NDraws',p.Results.FinalDraws);
groupDiagnostics=maxpoolfit.groupDiagnostics(cache,prediction);

result=struct('developmentFit',developmentFit,'finalFit',finalFit, ...
  'prediction',prediction,'groupDiagnostics',groupDiagnostics, ...
  'cachePath',string(cachePath),'createdAt',datetime('now'));
result.adaptationSelection=adaptationSelection(cfg,p.Results.AdaptationProfile);
if ~isempty(p.Results.AdaptationProfile)
  result.adaptationProfile=p.Results.AdaptationProfile;
end
if isfield(cache,'sourceManifest'), result.sourceManifest=cache.sourceManifest; end
if isfield(cache,'updatedAt'), result.cacheUpdatedAt=cache.updatedAt; end
if isfield(cache,'extractorVersion'), result.extractorVersion=cache.extractorVersion; end
if cfg.adaptationGainStrength>1e-12 && cfg.decisionRule=="absolute-step-max"
  resultStem='MaxPoolingIncrementAdaptationGainFit';
  noiseFile='MaxPoolingIncrementAdaptationGainNoiseSDByDomain.pdf';
elseif cfg.adaptationGainStrength>1e-12
  resultStem='MaxPoolingDirectionalIncrementAdaptationGainFit';
  noiseFile='MaxPoolingDirectionalIncrementAdaptationGainNoiseSDByDomain.pdf';
elseif cfg.adaptationExponent>1e-12 && cfg.decisionRule=="absolute-step-max"
  resultStem='MaxPoolingPowerAdaptationFit';
  noiseFile='MaxPoolingPowerAdaptationNoiseSDByDomain.pdf';
elseif cfg.adaptationExponent>1e-12
  resultStem='MaxPoolingDirectionalChangePowerAdaptationFit';
  noiseFile='MaxPoolingDirectionalChangePowerAdaptationNoiseSDByDomain.pdf';
elseif cfg.decisionRule=="directional-change-max"
  resultStem='MaxPoolingDirectionalChangeFit';
  noiseFile='MaxPoolingDirectionalChangeNoiseSDByDomain.pdf';
else
  resultStem='MaxPoolingFit';
  noiseFile='MaxPoolingNoiseSDByDomain.pdf';
end
save(fullfile(outputFolder,[resultStem 'Results.mat']),'result','-v7.3');

fig=figure('Color','w');
domainOrder=categorical(finalFit.domainFit.domain,finalFit.domainFit.domain,'Ordinal',true);
b=bar(domainOrder,[developmentFit.domainFit.noiseSDPC finalFit.domainFit.noiseSDPC]);
b(2).FaceAlpha=0.75;
ylabel('Fitted internal-noise SD (% coherence)');
title(sprintf(['Hard-max %s: domain noise, shared adaptation ' ...
  'a = %.3f, q = %.3f, b = %.3f'],strrep(char(cfg.decisionRule),'-',' '), ...
  finalFit.adaptation,cfg.adaptationExponent,cfg.adaptationGainStrength));
grid on; box off; legend({'Development','Final'},'Location','best');
exportgraphics(fig,fullfile(outputFolder,noiseFile), ...
  'ContentType','vector');
close(fig);
end

function selection=adaptationSelection(cfg,profileInput)
if cfg.decisionRule=="directional-change-max" && cfg.adaptationGainStrength>0
  equation="-a*Apre + exp(-b*Ppre/P0)*(Rstep-Ppre)";
elseif cfg.adaptationGainStrength>0
  equation="Ppre-a*Apre + exp(-b*Ppre/P0)*(Rstep-Ppre)";
elseif cfg.decisionRule=="directional-change-max" && cfg.adaptationExponent>0
  equation="Rstep - Rpre - a*Ppre*(Ppre/P0)^q";
elseif cfg.decisionRule=="directional-change-max"
  equation="Rstep - Rpre - a*Ppre";
elseif cfg.adaptationExponent>0
  equation="Rstep - a*Ppre*(Ppre/P0)^q";
else
  equation="Rstep - a*Ppre";
end
selection=struct('selectedAdaptation',cfg.adaptation, ...
  'adaptationExponent',cfg.adaptationExponent, ...
  'adaptationGainStrength',cfg.adaptationGainStrength, ...
  'adaptationReferencePC',cfg.adaptationReferencePC, ...
  'modelEquation',equation,'decisionRule',cfg.decisionRule, ...
  'sharedAcrossDomains',true,'sharedAcrossAnimals',true, ...
  'source',string(cfg.adaptationSelection));
if isempty(profileInput), return; end
assert(isfield(profileInput,'profile') && istable(profileInput.profile) && ...
  all(ismember({'adaptation','logLikelihood'}, ...
  profileInput.profile.Properties.VariableNames)), ...
  'AdaptationProfile must contain a profile table with adaptation and logLikelihood.');
P=profileInput.profile;
[bestLogLikelihood,bestRow]=max(P.logLikelihood);
selection.source="attached profile likelihood";
selection.profileBestGridAdaptation=P.adaptation(bestRow);
selection.profileBestLogLikelihood=bestLogLikelihood;
if isfield(profileInput,'nNoiseDraws')
  selection.profileNoiseDraws=profileInput.nNoiseDraws;
end
selection.profileGrid=[min(P.adaptation) max(P.adaptation)];
end
