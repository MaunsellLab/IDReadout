function result = profileNoiseCorrelationLength(cache,baseFit,varargin)
%PROFILENOISECORRELATIONLENGTH Sweep directional-noise correlation length.
%
% Group adaptations remain fixed. One noise SD per design group is
% reoptimized at every correlation length. Changed-side gain errors are
% unweighted across domain-animal-offset cells and exclude cells supported
% by fewer than MinimumGainSessions sessions.

p=inputParser;
p.addParameter('CorrelationLengthGridDeg',[5 10 15 20 30 45 60 90], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>0));
p.addParameter('NDraws',64, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('MinimumGainSessions',3, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.addParameter('ComputeGainDiagnostics',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
required={'trialGroup','trialMask','adaptationByGroup','config'};
assert(all(isfield(baseFit,required)), ...
  'A completed design-group fit is required.');
assert(numel(baseFit.trialGroup)==height(cache.trials), ...
  'Design-group labels do not match the cache.');

grid=unique(double(p.Results.CorrelationLengthGridDeg(:)),'stable');
trialGroup=string(baseFit.trialGroup(:));
adaptationByGroup=baseFit.adaptationByGroup;
workCache=cache; workCache.trials.domain=trialGroup;
nGrid=numel(grid); fits=cell(nGrid,1);
overallRows=cell(nGrid,8); groupRows=cell(0,16); gainRows=table();

for k=1:nGrid
  cfg=baseFit.config;
  cfg.noiseCorrLengthDeg=grid(k);
  cfg.modelVersion=string(baseFit.modelVersion)+"-corr-length-profile";
  cfg.verbose=false;
  fprintf('Correlation length %g deg (%d/%d)\n',grid(k),k,nGrid);
  fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
    adaptationByGroup,'TrialMask',baseFit.trialMask, ...
    'NDraws',p.Results.NDraws);
  fits{k}=fit;

  if p.Results.ComputeGainDiagnostics
    prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg, ...
      'NDraws',p.Results.NDraws);
    comparison=maxpoolfit.gainComparison(workCache,prediction,prediction, ...
      'TrialMask',baseFit.trialMask);
    G=comparison.summary;
    G=G(ismember(string(G.source), ...
      ["observed" "primary domain noise"]),:);
    G.correlationLengthDeg=repmat(grid(k),height(G),1);
    Gout=G;
    Gout.source(string(Gout.source)=="primary domain noise")= ...
      "primary design-group noise";
    gainRows=[gainRows;Gout]; %#ok<AGROW>
  else
    G=table();
  end

  groupMetric=nan(height(fit.groupFit),7);
  for g=1:height(fit.groupFit)
    label=fit.groupFit.analysisGroup(g);
    if p.Results.ComputeGainDiagnostics
      [metric,nConditions]=gainMetrics(G,label,p.Results.MinimumGainSessions);
      groupMetric(g,:)=[nConditions metric];
    end
    groupRows(end+1,:)={grid(k),label,fit.groupFit.domain(g), ...
      fit.groupFit.nTrials(g),fit.groupFit.nSessions(g), ...
      fit.groupFit.adaptation(g),fit.groupFit.noiseSDPC(g), ...
      fit.groupFit.logLikelihood(g),fit.groupFit.evalCount(g), ...
      groupMetric(g,1),groupMetric(g,2),groupMetric(g,3), ...
      groupMetric(g,4),groupMetric(g,5),groupMetric(g,6), ...
      groupMetric(g,7)}; %#ok<AGROW>
  end
  overallMetric=aggregateMetrics(groupMetric);
  overallRows(k,:)={grid(k),fit.logLikelihood,sum(fit.groupFit.evalCount), ...
    sum(groupMetric(:,1),'omitnan'),overallMetric(1),overallMetric(2), ...
    overallMetric(3),overallMetric(4)};
end

groupProfile=cell2table(groupRows,'VariableNames',{ ...
  'correlationLengthDeg','analysisGroup','domain','nTrials','nSessions', ...
  'adaptation','noiseSDPC','logLikelihood','evalCount','nGainConditions', ...
  'driftGainBias','probeGainBias','driftGainRMSE','probeGainRMSE', ...
  'absoluteGainRMSE','ratioRMSE'});
groupProfile.analysisGroup=string(groupProfile.analysisGroup);
groupProfile.domain=string(groupProfile.domain);
overallProfile=cell2table(overallRows,'VariableNames',{ ...
  'correlationLengthDeg','logLikelihood','evalCount','nGainConditions', ...
  'driftGainRMSE','probeGainRMSE','absoluteGainRMSE','ratioRMSE'});

overallProfile.deltaLogLFromBest=overallProfile.logLikelihood- ...
  max(overallProfile.logLikelihood);
groupProfile.deltaLogLFromGroupBest=nan(height(groupProfile),1);
groups=unique(groupProfile.analysisGroup,'stable');
for g=1:numel(groups)
  q=groupProfile.analysisGroup==groups(g);
  groupProfile.deltaLogLFromGroupBest(q)=groupProfile.logLikelihood(q)- ...
    max(groupProfile.logLikelihood(q));
end

result=struct('overallProfile',overallProfile,'groupProfile',groupProfile, ...
  'gainProfile',gainRows,'fits',{fits}, ...
  'correlationLengthGridDeg',grid,'nNoiseDraws',p.Results.NDraws, ...
  'minimumGainSessions',p.Results.MinimumGainSessions, ...
  'gainErrorWeighting',"unweighted across retained domain-animal-offset cells", ...
  'fixedAdaptationByGroup',adaptationByGroup, ...
  'baseModelVersion',baseFit.modelVersion,'createdAt',datetime('now'));
end

function [metric,nConditions]=gainMetrics(G,label,minimumSessions)
O=sortrows(G(string(G.domain)==label & string(G.source)=="observed",:), ...
  {'animal','probeOffsetDeg'});
M=sortrows(G(string(G.domain)==label & ...
  string(G.source)=="primary domain noise",:),{'animal','probeOffsetDeg'});
assert(height(O)==height(M) && ...
  isequal(string(O.animal),string(M.animal)) && ...
  isequal(O.probeOffsetDeg,M.probeOffsetDeg), ...
  'Observed and predicted gain cells do not align for %s.',label);
q=O.nSessions>=minimumSessions;
nConditions=sum(q);
if ~any(q), metric=nan(1,6); return; end
driftResidual=M.changedDriftGain(q)-O.changedDriftGain(q);
probeResidual=M.changedProbeGain(q)-O.changedProbeGain(q);
ratioResidual=M.changedProbeOverDrift(q)-O.changedProbeOverDrift(q);
metric=[mean(driftResidual) mean(probeResidual) ...
  sqrt(mean(driftResidual.^2)) sqrt(mean(probeResidual.^2)) ...
  sqrt(mean([driftResidual.^2;probeResidual.^2])) ...
  sqrt(mean(ratioResidual.^2))];
end

function metric=aggregateMetrics(groupMetric)
% Group rows contain counts followed by bias, component RMSEs, joint RMSE,
% and ratio RMSE. Reconstruct pooled sums without weighting by trial count.
q=isfinite(groupMetric(:,1)) & groupMetric(:,1)>0;
if ~any(q), metric=nan(1,4); return; end
n=groupMetric(q,1);
drift=sqrt(sum(n.*groupMetric(q,4).^2)/sum(n));
probe=sqrt(sum(n.*groupMetric(q,5).^2)/sum(n));
absolute=sqrt((sum(n.*groupMetric(q,4).^2)+ ...
  sum(n.*groupMetric(q,5).^2))/(2*sum(n)));
ratio=sqrt(sum(n.*groupMetric(q,7).^2)/sum(n));
metric=[drift probe absolute ratio];
end
