function performanceMatched = runPerformanceMatchedMaxPoolingFit( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNPERFORMANCEMATCHEDMAXPOOLINGFIT Fit condition noise to observed accuracy.
%
% performanceMatched = runPerformanceMatchedMaxPoolingFit( ...
%     cachePath,result,outputFolder)

p=inputParser;
p.addParameter('NDraws',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=2));
p.addParameter('GroupVariables', ...
  {'domain','animal','probeOffsetDeg','stepMagnitudePC'},@iscellstr);
p.addParameter('CalibrationGridSize',17, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=5&&mod(x,1)==0);
p.addParameter('MatchTolerance',0.005, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>0);
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
  'A completed primary fit is required.');
isGrouped=isfield(result.finalFit,'adaptationByGroup');
assert((isGrouped || isfinite(result.finalFit.adaptation)), ...
  'Primary fit has no valid adaptation specification.');
assert(isfield(result.finalFit,'config') && ...
  string(result.finalFit.modelVersion)==string(result.finalFit.config.modelVersion), ...
  'The supplied fit has inconsistent model configuration metadata.');
cfg=result.finalFit.config;
workCache=cache;
adaptationByTrial=[];
groupVariables=p.Results.GroupVariables;
if isGrouped
  assert(numel(result.finalFit.trialGroup)==height(cache.trials), ...
    'Grouped fit labels do not match the cache.');
  workCache.trials.domain=result.finalFit.trialGroup;
  adaptationByTrial=result.finalFit.trialAdaptation;
end
if isempty(p.Results.NDraws)
  nDraws=result.finalFit.nNoiseDraws;
else
  nDraws=p.Results.NDraws;
end

[fit,prediction]=maxpoolfit.fitPerformanceMatched(workCache,cfg, ...
  'TrialMask',result.finalFit.trialMask, ...
  'AdaptationByTrial',adaptationByTrial, ...
  'NDraws',nDraws, ...
  'GroupVariables',groupVariables, ...
  'CalibrationGridSize',p.Results.CalibrationGridSize, ...
  'MatchTolerance',p.Results.MatchTolerance);
calibration=maxpoolfit.calibrationDiagnostics(workCache,prediction);
if isGrouped
  fit.adaptationByGroup=result.finalFit.adaptationByGroup;
  prediction.adaptationByGroup=result.finalFit.adaptationByGroup;
end

% Put the original domain-noise prediction on the same condition rows.
idx=prediction.trialIdx;
T=workCache.trials(idx,:);
[groupID,keys]=findgroups(T(:,groupVariables));
assert(isequal(keys,fit.groupFit(:,groupVariables)), ...
  'Internal condition grouping mismatch.');
primaryP=result.prediction.pRDK1(idx);
assert(all(isfinite(primaryP)),'Primary predictions are incomplete.');
change1=double(T.changeRDK==1);
primaryCorrect=change1.*primaryP+(1-change1).*(1-primaryP);
primaryPredictedAccuracy=splitapply(@mean,primaryCorrect,groupID);
fit.groupFit.primaryPredictedAccuracy=primaryPredictedAccuracy;
fit.groupFit.primaryAccuracyResidual=fit.groupFit.observedAccuracy- ...
  primaryPredictedAccuracy;

performanceMatched=struct('fit',fit,'prediction',prediction, ...
  'calibration',calibration,'primaryModelVersion',result.finalFit.modelVersion, ...
  'cachePath',string(cachePath),'createdAt',datetime('now'));
save(fullfile(outputFolder,'PerformanceMatchedMaxPoolingFit.mat'), ...
  'performanceMatched','-v7.3');
writetable(fit.groupFit, ...
  fullfile(outputFolder,'PerformanceMatchedMaxPoolingFit_Conditions.csv'));
writetable(calibration.groupSummary, ...
  fullfile(outputFolder,'PerformanceMatchedMaxPoolingFit_Groups.csv'));
plotAccuracyResiduals(fit.groupFit,outputFolder);
end

function plotAccuracyResiduals(G,outputFolder)
labels=string(G.domain)+" - "+string(G.animal)+" - "+ ...
  string(G.probeOffsetDeg)+" deg - "+string(G.stepMagnitudePC)+"%";
x=categorical(labels,labels,'Ordinal',true);
fig=figure('Color','w','Position',[100 100 900 900]);
ax=axes(fig);
barh(ax,x,100.*[G.primaryAccuracyResidual G.accuracyResidual]);
xline(0,'k--');
xlabel('Observed - predicted accuracy (percentage points)');
ylabel('Condition'); grid on; box off;
ax.FontSize=7;
legend({'Primary domain-noise model','Performance matched'},'Location','best');
title('Accuracy residuals before and after condition calibration');
exportgraphics(fig,fullfile(outputFolder, ...
  'PerformanceMatchedMaxPoolingFit_AccuracyResiduals.pdf'),'ContentType','vector');
close(fig);
end
