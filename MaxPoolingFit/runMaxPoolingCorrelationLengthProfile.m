function correlationProfile = runMaxPoolingCorrelationLengthProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGCORRELATIONLENGTHPROFILE Profile directional-noise correlation.

p=inputParser;
p.addParameter('CorrelationLengthGridDeg',[5 10 15 20 30 45 60 90], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>0));
p.addParameter('NDraws',64, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('MinimumGainSessions',3, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
result=loadInput(resultInput,'result');
assert(isfield(result,'finalFit')&& ...
  isfield(result.finalFit,'adaptationByGroup'), ...
  'A completed design-group fit is required.');

correlationProfile=maxpoolfit.profileNoiseCorrelationLength( ...
  cache,result.finalFit, ...
  'CorrelationLengthGridDeg',p.Results.CorrelationLengthGridDeg, ...
  'NDraws',p.Results.NDraws, ...
  'MinimumGainSessions',p.Results.MinimumGainSessions);
correlationProfile.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingCorrelationLengthProfile.mat'), ...
  'correlationProfile','-v7.3');
writetable(correlationProfile.overallProfile,fullfile(outputFolder, ...
  'MaxPoolingCorrelationLengthProfile_Overall.csv'));
writetable(correlationProfile.groupProfile,fullfile(outputFolder, ...
  'MaxPoolingCorrelationLengthProfile_Groups.csv'));
writetable(correlationProfile.gainProfile,fullfile(outputFolder, ...
  'MaxPoolingCorrelationLengthProfile_Gains.csv'));
plotProfile(correlationProfile,outputFolder);
end

function value=loadInput(input,name)
if isstruct(input), value=input; return; end
S=load(input,name); assert(isfield(S,name), ...
  'MAT file does not contain variable %s.',name); value=S.(name);
end

function plotProfile(P,outputFolder)
O=P.overallProfile; G=P.groupProfile;
groups=unique(G.analysisGroup,'stable'); colors=lines(numel(groups));
fig=figure('Color','w','Position',[80 80 1100 780]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
plot(ax,O.correlationLengthDeg,O.deltaLogLFromBest,'ko-', ...
  'MarkerFaceColor','k','LineWidth',1.3);
xline(ax,30,'k:','Current','HandleVisibility','off');
xlabel(ax,'Noise correlation length (deg)'); ylabel(ax,'Delta log likelihood');
title(ax,'Combined choice likelihood');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
for g=1:numel(groups)
  q=G.analysisGroup==groups(g);
  plot(ax,G.correlationLengthDeg(q),G.deltaLogLFromGroupBest(q),'o-', ...
    'Color',colors(g,:),'MarkerFaceColor',colors(g,:), ...
    'LineWidth',1.2,'DisplayName',groups(g));
end
xline(ax,30,'k:','Current','HandleVisibility','off');
xlabel(ax,'Noise correlation length (deg)'); ylabel(ax,'Delta log likelihood');
title(ax,'Group likelihood'); legend(ax,'Location','best');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
for g=1:numel(groups)
  q=G.analysisGroup==groups(g);
  plot(ax,G.correlationLengthDeg(q),G.noiseSDPC(q),'o-', ...
    'Color',colors(g,:),'MarkerFaceColor',colors(g,:), ...
    'LineWidth',1.2,'DisplayName',groups(g));
end
xline(ax,30,'k:','Current','HandleVisibility','off');
xlabel(ax,'Noise correlation length (deg)');
ylabel(ax,'Refitted noise SD (% coherence)'); title(ax,'Noise tradeoff');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
for g=1:numel(groups)
  q=G.analysisGroup==groups(g);
  plot(ax,G.correlationLengthDeg(q),G.probeGainRMSE(q),'o-', ...
    'Color',colors(g,:),'MarkerFaceColor',colors(g,:), ...
    'LineWidth',1.2,'DisplayName',groups(g));
end
xline(ax,30,'k:','Current','HandleVisibility','off');
xlabel(ax,'Noise correlation length (deg)');
ylabel(ax,'Changed-probe gain RMSE'); title(ax,'Gain geometry (unweighted)');

for ax=findall(fig,'Type','axes').'
  ax.XScale='log'; ax.XTick=O.correlationLengthDeg;
  ax.XTickLabel=string(O.correlationLengthDeg);
end
title(tl,sprintf(['Directional-noise correlation profile; %d draws; ' ...
  'gain cells require >=%d sessions'],P.nNoiseDraws,P.minimumGainSessions));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingCorrelationLengthProfile.pdf'),'ContentType','vector');
close(fig);
end
