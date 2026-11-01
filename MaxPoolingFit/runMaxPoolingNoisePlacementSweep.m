function noisePlacement = runMaxPoolingNoisePlacementSweep( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGNOISEPLACEMENTSWEEP Sweep pre-max noise and solve post-WTA flips.

p=inputParser;
p.addParameter('Domain',"IDR",@(x)isstring(x)||ischar(x));
p.addParameter('Animal',"Meetz",@(x)isstring(x)||ischar(x));
p.addParameter('ProbeOffsetDeg',25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('StepMagnitudePC',23,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('NoiseSDGridPC',0:1:24, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('MinimumTrials',500, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('MinimumSessions',5, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end

cache=maxpoolfit.loadCache(cachePath);
if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit'),'A completed primary fit is required.');

noisePlacement=maxpoolfit.noisePlacementSweep(cache,result.finalFit, ...
  'Domain',p.Results.Domain,'Animal',p.Results.Animal, ...
  'ProbeOffsetDeg',p.Results.ProbeOffsetDeg, ...
  'StepMagnitudePC',p.Results.StepMagnitudePC, ...
  'NoiseSDGridPC',p.Results.NoiseSDGridPC,'NDraws',p.Results.NDraws, ...
  'MinimumTrials',p.Results.MinimumTrials, ...
  'MinimumSessions',p.Results.MinimumSessions);
noisePlacement.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingNoisePlacementSweep.mat'), ...
  'noisePlacement','-v7.3');
writetable(noisePlacement.sweep,fullfile(outputFolder, ...
  'MaxPoolingNoisePlacementSweep.csv'));
writetable(noisePlacement.targets,fullfile(outputFolder, ...
  'MaxPoolingNoisePlacementTargets.csv'));
plotOverview(noisePlacement,outputFolder);
plotAbsoluteGains(noisePlacement,outputFolder);
end

function plotOverview(D,outputFolder)
S=D.sweep; targetDrift=D.targets.unchangedOverChanged(D.targets.component=="drift");
targetProbe=D.targets.unchangedOverChanged(D.targets.component=="probe");
fig=figure('Color','w','Position',[100 100 1050 780]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
plot(ax,S.noiseSDPC,S.upstreamDriftUnchangedOverChanged,'o-', ...
  'LineWidth',1.3,'DisplayName','Drift');
plot(ax,S.noiseSDPC,S.upstreamProbeUnchangedOverChanged,'s-', ...
  'LineWidth',1.3,'DisplayName','Probe');
yline(ax,targetDrift,'--','Observed drift','HandleVisibility','off');
yline(ax,targetProbe,':','Observed probe','HandleVisibility','off');
markNoise(ax,D); xlabel(ax,'Pre-max noise SD (% coherence)');
ylabel(ax,'Unchanged / changed gain'); title(ax,'Upstream gain ratios');
legend(ax,'Location','best');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
plot(ax,S.noiseSDPC,S.upstreamAccuracy,'o-','LineWidth',1.3, ...
  'DisplayName','Before downstream flip');
plot(ax,S.noiseSDPC,S.postFlipAccuracy,'s--','LineWidth',1.3, ...
  'DisplayName','After fitted flip');
yline(ax,D.observedAccuracy,'k:','Observed','HandleVisibility','off');
markNoise(ax,D); xlabel(ax,'Pre-max noise SD (% coherence)');
ylabel(ax,'Accuracy'); ylim(ax,[0.45 1]); title(ax,'Accuracy');
legend(ax,'Location','best');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
plot(ax,S.noiseSDPC,100.*S.requiredFlipRate,'o-','LineWidth',1.3);
markNoise(ax,D); xlabel(ax,'Pre-max noise SD (% coherence)');
ylabel(ax,'Required post-selection flips (%)'); title(ax,'Accuracy-matching flip rate');

ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
plot(ax,S.noiseSDPC,S.upstreamLogLoss,'o-','LineWidth',1.3, ...
  'DisplayName','Before downstream flip');
plot(ax,S.noiseSDPC,S.postFlipLogLoss,'s--','LineWidth',1.3, ...
  'DisplayName','After fitted flip');
markNoise(ax,D); xlabel(ax,'Pre-max noise SD (% coherence)');
ylabel(ax,'Log loss per trial'); title(ax,'Trial-level prediction');
legend(ax,'Location','best');

title(tl,selectionTitle(D,"Noise placement: ratio first, accuracy second"));
exportgraphics(fig,fullfile(outputFolder,'MaxPoolingNoisePlacementSweep.pdf'), ...
  'ContentType','vector'); close(fig);
end

function plotAbsoluteGains(D,outputFolder)
S=D.sweep; fig=figure('Color','w','Position',[100 100 1050 470]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
components=["Drift" "Probe"];
for k=1:2
  ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
  changed="postFlipChanged"+components(k)+"CoefficientPerPC";
  unchanged="postFlipUnchanged"+components(k)+"CoefficientPerPC";
  q=D.targets.component==lower(components(k));
  plot(ax,S.noiseSDPC,S.(changed),'o-','LineWidth',1.3, ...
    'DisplayName','Changed RDK');
  plot(ax,S.noiseSDPC,S.(unchanged),'s--','LineWidth',1.3, ...
    'DisplayName','-Unchanged RDK');
  yline(ax,D.targets.changedCoefficientPerPC(q),'--','Observed changed', ...
    'HandleVisibility','off');
  yline(ax,D.targets.unchangedCoefficientPerPC(q),':','Observed -unchanged', ...
    'HandleVisibility','off');
  markNoise(ax,D); xlabel(ax,'Pre-max noise SD (% coherence)');
  ylabel(ax,'Choice coefficient per % coherence'); title(ax,components(k));
  if k==1, legend(ax,'Location','best'); end
end
title(tl,selectionTitle(D,"Absolute gains after accuracy-matching flips"));
exportgraphics(fig,fullfile(outputFolder,'MaxPoolingNoisePlacementGains.pdf'), ...
  'ContentType','vector'); close(fig);
end

function markNoise(ax,D)
xline(ax,D.currentNoiseSDPC,'k-.','Current fit','HandleVisibility','off');
if ~isempty(D.bestFeasibleRatioRow)
  xline(ax,D.bestFeasibleRatioRow.noiseSDPC,'Color',[0.1 0.55 0.2], ...
    'LineStyle','--','Label','Best feasible ratio','HandleVisibility','off');
end
end

function value=selectionTitle(D,prefix)
Q=D.selection;
value=sprintf('%s\n%s - %s; %g deg probe; %g%% step; n=%d, sessions=%d', ...
  prefix,Q.domain,Q.animal,Q.probeOffsetDeg,Q.stepMagnitudePC,Q.nTrials,Q.nSessions);
end
