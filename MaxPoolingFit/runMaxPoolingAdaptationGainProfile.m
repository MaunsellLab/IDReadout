function gainProfile = runMaxPoolingAdaptationGainProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGADAPTATIONGAINPROFILE Test exponential adaptation of response gain.

p=inputParser;
p.addParameter('GainGrid',0:0.1:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('Adaptation',[], ...
  @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1));
p.addParameter('ReferencePC',40, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit')&&isfinite(result.finalFit.adaptation), ...
  'A completed primary fit is required.');
cache=maxpoolfit.loadCache(cachePath);
cfg=result.finalFit.config;
if isempty(p.Results.Adaptation), adaptation=result.finalFit.adaptation;
else, adaptation=double(p.Results.Adaptation); end
cfg.adaptation=adaptation;
cfg=maxpoolfit.configureAdaptationExponent(cfg,0,p.Results.ReferencePC);
cfg=maxpoolfit.configureAdaptationGain(cfg,0,p.Results.ReferencePC);
gainProfile=maxpoolfit.profileAdaptationGain(cache,cfg, ...
  'TrialMask',result.finalFit.trialMask,'NDraws',p.Results.NDraws, ...
  'GainGrid',p.Results.GainGrid,'Adaptation',adaptation, ...
  'ReferencePC',p.Results.ReferencePC);
stem='MaxPoolingIncrementAdaptationGainProfile';
save(fullfile(outputFolder,[stem '.mat']),'gainProfile','-v7.3');

P=gainProfile.profile;
exportTable=P(:,{'adaptation','adaptationGainStrength','logLikelihood', ...
  'deltaLogLFromBest','deltaLogLFromLinear'});
for id=1:numel(gainProfile.domains)
  suffix=matlab.lang.makeValidName(gainProfile.domains(id));
  exportTable.("noiseSDPC_"+suffix)=P.noiseSDPC(:,id);
  exportTable.("logLikelihood_"+suffix)=P.domainLogLikelihood(:,id);
  exportTable.("evalCount_"+suffix)=P.evalCount(:,id);
end
writetable(exportTable,fullfile(outputFolder,[stem '.csv']));

fig=figure('Color','w','Position',[100 100 900 450]);
plot(P.adaptationGainStrength,P.deltaLogLFromBest,'ko-', ...
  'LineWidth',1.5,'MarkerFaceColor','k','DisplayName','All domains');
hold on; colors=lines(numel(gainProfile.domains));
for id=1:numel(gainProfile.domains)
  y=P.domainLogLikelihood(:,id);
  plot(P.adaptationGainStrength,y-max(y),'o-','Color',colors(id,:), ...
    'LineWidth',1.2,'DisplayName',gainProfile.domains(id));
end
yline(0,'k:','HandleVisibility','off'); grid on; box off;
xlabel('Adaptation gain strength, b');
ylabel('\Delta log likelihood from profile maximum');
title(sprintf(['Increment-only prestep-dependent gain: a = %.3f, ' ...
  'P_0 = %.1f%%'],adaptation,p.Results.ReferencePC));
legend('Location','best');
exportgraphics(fig,fullfile(outputFolder,[stem '.pdf']),'ContentType','vector');
close(fig);
end
