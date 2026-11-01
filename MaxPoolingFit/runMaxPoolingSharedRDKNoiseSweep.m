function sharedNoise = runMaxPoolingSharedRDKNoiseSweep( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGSHAREDRDKNOISESWEEP Test noise shared between the two RDKs.

if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath); result=loadResult(resultInput);
assert(isfield(result,'finalFit'),'A completed primary fit is required.');
sharedNoise=maxpoolfit.sharedRDKNoiseSweep(cache,result.finalFit,varargin{:});
sharedNoise.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingSharedRDKNoiseSweep.mat'), ...
  'sharedNoise','-v7.3');
writetable(sharedNoise.sweep,fullfile(outputFolder, ...
  'MaxPoolingSharedRDKNoiseSweep.csv'));
writetable(sharedNoise.targets,fullfile(outputFolder, ...
  'MaxPoolingSharedRDKNoiseTargets.csv'));
writetable(sharedNoise.bestFeasibleByCorrelation,fullfile(outputFolder, ...
  'MaxPoolingSharedRDKNoiseBestByCorrelation.csv'));
plotSweep(sharedNoise,outputFolder);
end

function result=loadResult(input)
if isstruct(input), result=input; return; end
S=load(input,'result');
assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
result=S.result;
end

function plotSweep(D,outputFolder)
S=D.sweep; rhos=unique(S.rdkNoiseCorrelation,'stable'); colors=parula(numel(rhos));
targetDrift=D.targets.unchangedOverChanged(D.targets.component=="drift");
targetProbe=D.targets.unchangedOverChanged(D.targets.component=="probe");
fig=figure('Color','w','Position',[100 100 1100 800]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); plotLines(ax,S,rhos,colors, ...
  'postFlipDriftUnchangedOverChanged');
yline(ax,targetDrift,'k--','Observed','HandleVisibility','off');
ylabel(ax,'Unchanged / changed drift gain'); title(ax,'Drift ratio after flips');
ax=nexttile(tl); plotLines(ax,S,rhos,colors, ...
  'postFlipProbeUnchangedOverChanged');
yline(ax,targetProbe,'k--','Observed','HandleVisibility','off');
ylabel(ax,'Unchanged / changed probe gain'); title(ax,'Probe ratio after flips');
ax=nexttile(tl); plotLines(ax,S,rhos,colors,'upstreamAccuracy');
yline(ax,D.observedAccuracy,'k--','Observed','HandleVisibility','off');
ylabel(ax,'Accuracy before downstream flips'); title(ax,'Upstream accuracy');
ax=nexttile(tl); plotLines(ax,S,rhos,colors,'requiredFlipRate');
ylabel(ax,'Required downstream flip rate'); title(ax,'Accuracy matching');
Q=D.selection;
title(tl,sprintf('Shared RDK noise: %s - %s; %g deg; %g%% step; paired draws=%d', ...
  Q.domain,Q.animal,Q.probeOffsetDeg,Q.stepMagnitudePC,D.nNoiseDraws));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingSharedRDKNoiseSweep.pdf'),'ContentType','vector'); close(fig);
end

function plotLines(ax,S,rhos,colors,variable)
hold(ax,'on'); grid(ax,'on'); box(ax,'off');
for k=1:numel(rhos)
  q=S.rdkNoiseCorrelation==rhos(k); R=sortrows(S(q,:),'noiseSDPC');
  plot(ax,R.noiseSDPC,R.(variable),'o-','Color',colors(k,:), ...
    'LineWidth',1.2,'DisplayName',sprintf('rho = %.2g',rhos(k)));
end
xlabel(ax,'Pre-max noise SD (% coherence)'); legend(ax,'Location','best');
end
