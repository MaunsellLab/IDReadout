function adaptationProfile = runMaxPoolingAdaptationProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGADAPTATIONPROFILE Profile Rstep-a*Ppre with domain noise.

p=inputParser;
p.addParameter('AdaptationGrid',0:0.1:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('DecisionRule',"", ...
  @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit') && isfinite(result.finalFit.adaptation), ...
  'A completed full-prestep primary fit is required.');
cache=maxpoolfit.loadCache(cachePath);
cfg=result.finalFit.config;
if ~isfield(cfg,'decisionRule'), cfg.decisionRule="absolute-step-max"; end
if strlength(string(p.Results.DecisionRule))>0
  cfg=maxpoolfit.configureDecisionRule(cfg,p.Results.DecisionRule);
else
  cfg=maxpoolfit.configureDecisionRule(cfg,cfg.decisionRule);
end
adaptationProfile=maxpoolfit.profileAdaptation(cache,cfg, ...
  'TrialMask',result.finalFit.trialMask, ...
  'NDraws',p.Results.NDraws, ...
  'AdaptationGrid',p.Results.AdaptationGrid);
if cfg.decisionRule=="directional-change-max"
  stem='MaxPoolingDirectionalChangeAdaptationProfile';
else
  stem='MaxPoolingAdaptationProfile';
end
save(fullfile(outputFolder,[stem '.mat']), ...
  'adaptationProfile','-v7.3');

P=adaptationProfile.profile;
exportTable=P(:,{'adaptation','logLikelihood','deltaLogLFromBest', ...
  'deltaLogLFromZero'});
for id=1:numel(adaptationProfile.domains)
  suffix=matlab.lang.makeValidName(adaptationProfile.domains(id));
  exportTable.("noiseSDPC_"+suffix)=P.noiseSDPC(:,id);
  exportTable.("logLikelihood_"+suffix)=P.domainLogLikelihood(:,id);
end
writetable(exportTable,fullfile(outputFolder,[stem '.csv']));

fig=figure('Color','w','Position',[100 100 900 450]);
plot(P.adaptation,P.deltaLogLFromBest,'ko-', ...
  'LineWidth',1.5,'MarkerFaceColor','k','DisplayName','All domains');
hold on;
colors=lines(numel(adaptationProfile.domains));
for id=1:numel(adaptationProfile.domains)
  y=P.domainLogLikelihood(:,id);
  plot(P.adaptation,y-max(y),'o-','Color',colors(id,:), ...
    'LineWidth',1.2,'DisplayName',adaptationProfile.domains(id));
end
yline(0,'k:','HandleVisibility','off'); grid on; box off;
xlabel('Adaptation, a'); ylabel('\Delta log likelihood from profile maximum');
if cfg.decisionRule=="directional-change-max"
  profileTitle='Directional-change profile: R_{step} - R_{pre} - aP_{pre}';
else
  profileTitle='Full-prestep adaptation profile: R_{step} - aP_{pre}';
end
title(profileTitle);
legend('Location','best');
exportgraphics(fig,fullfile(outputFolder,[stem '.pdf']),'ContentType','vector');
close(fig);
end
