function guessProfile = runMaxPoolingGuessCutoffProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGGUESSCUTOFFPROFILE Profile a low-confidence guess criterion.

p=inputParser;
p.addParameter('KappaGrid',[-Inf 0:0.25:3], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('GuessP_RDK1',0.5, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end

if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit') && isfield(result.finalFit,'config'), ...
  'A completed absolute-step-max fit is required.');
cfg=result.finalFit.config;
if ~isfield(cfg,'decisionRule'), cfg.decisionRule="absolute-step-max"; end
assert(string(cfg.decisionRule)=="absolute-step-max", ...
  'The guess-cutoff diagnostic must start from the absolute-step-max fit.');
cache=maxpoolfit.loadCache(cachePath);
guessProfile=maxpoolfit.profileGuessCutoff(cache,cfg, ...
  'TrialMask',result.finalFit.trialMask, ...
  'NDraws',p.Results.NDraws,'KappaGrid',p.Results.KappaGrid, ...
  'GuessP_RDK1',p.Results.GuessP_RDK1);

stem='MaxPoolingGuessCutoffProfile';
save(fullfile(outputFolder,[stem '.mat']),'guessProfile','-v7.3');
P=guessProfile.profile;
exportTable=P(:,{'kappa','logLikelihood','deltaLogLFromBest', ...
  'deltaLogLFromNoGuess'});
for id=1:numel(guessProfile.domains)
  suffix=matlab.lang.makeValidName(guessProfile.domains(id));
  exportTable.("noiseSDPC_"+suffix)=P.noiseSDPC(:,id);
  exportTable.("logLikelihood_"+suffix)=P.domainLogLikelihood(:,id);
  exportTable.("guessFraction_"+suffix)=P.guessFraction(:,id);
  exportTable.("predictedAccuracy_"+suffix)=P.predictedAccuracy(:,id);
end
writetable(exportTable,fullfile(outputFolder,[stem '.csv']));
writetable(guessProfile.conditionSummary,fullfile(outputFolder, ...
  'MaxPoolingGuessCutoffConditions.csv'));
plotProfile(guessProfile,outputFolder);
end

function plotProfile(G,outputFolder)
P=G.profile; finite=isfinite(P.kappa);
fig=figure('Color','w','Position',[100 100 950 430]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);
plot(ax,P.kappa(finite),P.deltaLogLFromBest(finite),'ko-', ...
  'LineWidth',1.4,'MarkerFaceColor','k'); hold(ax,'on');
reference=P.deltaLogLFromBest(~finite);
if ~isempty(reference)
  yline(ax,reference(1),'k--','No guess','LabelHorizontalAlignment','left');
end
yline(ax,0,'k:','HandleVisibility','off'); grid(ax,'on'); box(ax,'off');
xlabel(ax,'Detection cutoff, \kappa (domain-noise SD)');
ylabel(ax,'\Delta log likelihood from profile maximum');
title(ax,'Choice likelihood');

ax=nexttile(tl); colors=lines(numel(G.domains)); hold(ax,'on');
for id=1:numel(G.domains)
  plot(ax,P.kappa(finite),100.*P.guessFraction(finite,id),'o-', ...
    'Color',colors(id,:),'LineWidth',1.3,'DisplayName',G.domains(id));
end
grid(ax,'on'); box(ax,'off'); xlabel(ax,'Detection cutoff, \kappa');
ylabel(ax,'Mean guess probability (%)'); title(ax,'Guess-state occupancy');
legend(ax,'Location','best');
title(tl,sprintf('Detection-race guess profile; adaptation a = %.3f',G.adaptation));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingGuessCutoffProfile.pdf'),'ContentType','vector');
close(fig);
end
