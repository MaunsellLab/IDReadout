function calibration = runMaxPoolingCalibration(cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGCALIBRATION Compute, save, and plot fitted-trial calibration.
%
% calibration = runMaxPoolingCalibration(cachePath,result,outputFolder)
% calibration = runMaxPoolingCalibration(cachePath,resultPath,outputFolder)

p=inputParser;
p.addParameter('NBins',10,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
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
assert(isfield(result,'finalFit') && isfield(result,'prediction'), ...
  'Result must contain finalFit and prediction.');
isGrouped=isfield(result.finalFit,'adaptationByGroup');
assert((isGrouped || (isfield(result.finalFit,'domains') && ...
  isfinite(result.finalFit.adaptation))), ...
  'Calibration requires a completed shared or design-group fit.');
assert(isfield(result.finalFit,'config') && ...
  string(result.finalFit.modelVersion)==string(result.finalFit.config.modelVersion), ...
  'The supplied fit has inconsistent model configuration metadata.');

fitIdx=find(result.finalFit.trialMask);
assert(all(isfinite(result.prediction.pRDK1(fitIdx))), ...
  'Saved predictions do not cover every primary-fit trial.');
fitPrediction=result.prediction;
fitPrediction.trialIdx=fitIdx;
workCache=cache;
if isGrouped
  assert(numel(result.finalFit.trialGroup)==height(cache.trials), ...
    'Grouped fit labels do not match the cache.');
  workCache.trials.domain=result.finalFit.trialGroup;
end
calibration=maxpoolfit.calibrationDiagnostics(workCache,fitPrediction, ...
  'NBins',p.Results.NBins);
calibration.modelVersion=result.finalFit.modelVersion;
if isGrouped
  calibration.adaptationByGroup=result.finalFit.adaptationByGroup;
end
calibration.cachePath=string(cachePath);

save(fullfile(outputFolder,'MaxPoolingCalibration.mat'),'calibration','-v7.3');
writetable(calibration.groupSummary, ...
  fullfile(outputFolder,'MaxPoolingCalibration_Groups.csv'));
writetable(calibration.reliability, ...
  fullfile(outputFolder,'MaxPoolingCalibration_Reliability.csv'));
writetable(calibration.sessionSummary, ...
  fullfile(outputFolder,'MaxPoolingCalibration_Sessions.csv'));
writetable(calibration.conditionSummary, ...
  fullfile(outputFolder,'MaxPoolingCalibration_Conditions.csv'));

plotReliability(calibration,outputFolder);
plotGroupAccuracy(calibration,outputFolder);
plotSessionAccuracy(calibration,outputFolder);
end

function plotReliability(C,outputFolder)
G=C.groupSummary;
fig=figure('Color','w','Position',[100 100 1050 650]);
tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for g=1:height(G)
  ax=nexttile(tl);
  q=string(C.reliability.domain)==string(G.domain(g)) & ...
    string(C.reliability.animal)==string(G.animal(g));
  R=C.reliability(q,:);
  plot(ax,[0 1],[0 1],'k--','LineWidth',0.8); hold(ax,'on');
  errorbar(ax,R.meanPredictedPChooseRDK1,R.observedPChooseRDK1, ...
    R.observedPChooseRDK1-R.ciLow,R.ciHigh-R.observedPChooseRDK1, ...
    'o-','LineWidth',1.2,'MarkerFaceColor',[0.2 0.45 0.75]);
  axis(ax,'square'); xlim(ax,[0 1]); ylim(ax,[0 1]); grid(ax,'on'); box(ax,'off');
  title(ax,sprintf('%s - %s',G.domain(g),G.animal(g)));
  xlabel(ax,'Predicted P(choose RDK 1)'); ylabel(ax,'Observed P(choose RDK 1)');
end
title(tl,'Choice-probability calibration by domain and animal');
exportgraphics(fig,fullfile(outputFolder,'MaxPoolingCalibration_Reliability.pdf'), ...
  'ContentType','vector');
close(fig);
end

function plotGroupAccuracy(C,outputFolder)
G=C.groupSummary;
labels=string(G.domain)+" - "+string(G.animal);
x=categorical(labels,labels,'Ordinal',true);
fig=figure('Color','w','Position',[100 100 1050 450]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);
bar(ax,x,100.*[G.observedAccuracy G.predictedAccuracy]);
ylabel(ax,'Correct choices (%)'); grid(ax,'on'); box(ax,'off');
legend(ax,{'Observed','Predicted'},'Location','best');
title(ax,'Behavioral accuracy');
ax=nexttile(tl);
bar(ax,x,100.*G.accuracyResidual);
yline(ax,0,'k--'); ylabel(ax,'Observed - predicted (percentage points)');
grid(ax,'on'); box(ax,'off'); title(ax,'Accuracy residual');
title(tl,'Group-level calibration');
exportgraphics(fig,fullfile(outputFolder,'MaxPoolingCalibration_Groups.pdf'), ...
  'ContentType','vector');
close(fig);
end

function plotSessionAccuracy(C,outputFolder)
S=C.sessionSummary; G=C.groupSummary;
fig=figure('Color','w','Position',[100 100 650 600]); ax=axes(fig); hold(ax,'on');
colors=lines(height(G));
for g=1:height(G)
  q=string(S.domain)==string(G.domain(g)) & string(S.animal)==string(G.animal(g));
  scatter(ax,100.*S.predictedAccuracy(q),100.*S.observedAccuracy(q),22, ...
    colors(g,:),'filled','MarkerFaceAlpha',0.65, ...
    'DisplayName',sprintf('%s - %s',G.domain(g),G.animal(g)));
end
allAccuracy=100.*[S.predictedAccuracy;S.observedAccuracy];
lo=max(0,5*floor(min(allAccuracy)/5)-5);
hi=min(100,5*ceil(max(allAccuracy)/5)+5);
if lo>=hi, lo=max(0,hi-10); end
lims=[lo hi]; plot(ax,lims,lims,'k--','HandleVisibility','off');
xlim(ax,lims); ylim(ax,lims); axis(ax,'square'); grid(ax,'on'); box(ax,'off');
xlabel(ax,'Predicted correct choices (%)'); ylabel(ax,'Observed correct choices (%)');
title(ax,'Session-level accuracy calibration'); legend(ax,'Location','best');
exportgraphics(fig,fullfile(outputFolder,'MaxPoolingCalibration_Sessions.pdf'), ...
  'ContentType','vector');
close(fig);
end
