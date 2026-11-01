function [fit,prediction] = fitPerformanceMatched(cache,cfg,varargin)
%FITPERFORMANCEMATCHED Match observed accuracy within experimental conditions.
%
% This is a diagnostic calibration fit, not the primary likelihood fit. One
% noise SD is chosen per condition to minimize squared error between observed
% and predicted accuracy while all max-pooling mechanisms remain fixed.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('AdaptationByTrial',[],@(x)isempty(x)||isnumeric(x));
p.addParameter('NDraws',cfg.finalNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('GroupVariables', ...
  {'domain','animal','probeOffsetDeg','stepMagnitudePC'},@iscellstr);
p.addParameter('CalibrationGridSize',17, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=5&&mod(x,1)==0);
p.addParameter('MatchTolerance',0.005, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
if ~isfield(cfg,'adaptationExponent'), cfg.adaptationExponent=0; end
if ~isfield(cfg,'adaptationGainStrength'), cfg.adaptationGainStrength=0; end
if ~isfield(cfg,'adaptationReferencePC'), cfg.adaptationReferencePC=40; end
maxpoolfit.validateCache(cache);
if isempty(p.Results.AdaptationByTrial)
  assert(isfield(cfg,'adaptation') && isnumeric(cfg.adaptation) && ...
    isscalar(cfg.adaptation) && isfinite(cfg.adaptation) && ...
    cfg.adaptation>=0 && cfg.adaptation<=1, ...
    'Performance matching requires a valid shared adaptation value.');
  adaptationByTrial=repmat(double(cfg.adaptation),height(cache.trials),1);
else
  adaptationByTrial=double(p.Results.AdaptationByTrial(:));
  assert(numel(adaptationByTrial)==height(cache.trials) && ...
    all(isfinite(adaptationByTrial)) && ...
    all(adaptationByTrial>=0 & adaptationByTrial<=1), ...
    'AdaptationByTrial must contain one finite [0,1] value per cached trial.');
end
if isempty(p.Results.TrialMask)
  mask=cfg.primarySelector(cache.trials);
else
  mask=false(height(cache.trials),1);
  if islogical(p.Results.TrialMask)
    assert(numel(p.Results.TrialMask)==height(cache.trials), ...
      'Logical TrialMask must have one value per cached trial.');
    mask=p.Results.TrialMask(:);
  else
    mask(p.Results.TrialMask)=true;
  end
end
groupVars=p.Results.GroupVariables;
assert(all(ismember(groupVars,cache.trials.Properties.VariableNames)), ...
  'One or more requested performance-matching group variables are absent from cache.trials.');
selectedIdx=find(mask);
assert(~isempty(selectedIdx),'No trials selected.');
[groupID,groupKeys]=findgroups(cache.trials(selectedIdx,groupVars));
nGroups=height(groupKeys);

z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
logGrid=linspace(log(cfg.noiseBoundsPC(1)),log(cfg.noiseBoundsPC(2)), ...
  p.Results.CalibrationGridSize);
pRDK1=nan(height(cache.trials),1);
nTrials=zeros(nGroups,1);
observedAccuracy=nan(nGroups,1);
predictedAccuracy=nan(nGroups,1);
accuracyResidual=nan(nGroups,1);
noiseSDPC=nan(nGroups,1);
logLikelihood=nan(nGroups,1);
minGridAccuracy=nan(nGroups,1);
maxGridAccuracy=nan(nGroups,1);
evalCount=zeros(nGroups,1);
matchStatus=strings(nGroups,1);

for g=1:nGroups
  idx=selectedIdx(groupID==g);
  nTrials(g)=numel(idx);
  observedCorrect=double(cache.trials.choiceRDK(idx)==cache.trials.changeRDK(idx));
  observedAccuracy(g)=mean(observedCorrect);
  change1=double(cache.trials.changeRDK(idx)==1);
  choice1=double(cache.trials.choiceRDK(idx)==1);
  groupAdaptation=unique(adaptationByTrial(idx));
  assert(numel(groupAdaptation)==1, ...
    'A performance-matching condition contains multiple adaptation values.');
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  if cfg.verbose
    fprintf('Performance match %d/%d (%d trials)\n',g,nGroups,nTrials(g));
  end

  gridAccuracy=nan(size(logGrid));
  for k=1:numel(logGrid)
    gridAccuracy(k)=accuracyAtSigma(exp(logGrid(k)),B,A,Rpre,change1, ...
      groupAdaptation,cfg,z1,z2);
  end
  evalCount(g)=numel(logGrid);
  minGridAccuracy(g)=min(gridAccuracy);
  maxGridAccuracy(g)=max(gridAccuracy);
  [~,bestGrid]=min(abs(gridAccuracy-observedAccuracy(g)));
  lo=max(1,bestGrid-1); hi=min(numel(logGrid),bestGrid+1);
  candidates=logGrid(bestGrid);
  if lo<hi
    objective=@(logSigma) (accuracyAtSigma(exp(logSigma), ...
      B,A,Rpre,change1,groupAdaptation,cfg,z1,z2)-observedAccuracy(g)).^2;
    [x,~,~,out]=fminbnd(objective,logGrid(lo),logGrid(hi), ...
      optimset('TolX',2e-3,'Display','off'));
    candidates=[candidates x logGrid(lo) logGrid(hi)]; %#ok<AGROW>
    evalCount(g)=evalCount(g)+out.funcCount;
  end
  candidateAccuracy=arrayfun(@(x)accuracyAtSigma(exp(x), ...
    B,A,Rpre,change1,groupAdaptation,cfg,z1,z2),candidates);
  evalCount(g)=evalCount(g)+numel(candidates);
  [~,best]=min(abs(candidateAccuracy-observedAccuracy(g)));
  noiseSDPC(g)=exp(candidates(best));
  predictedAccuracy(g)=candidateAccuracy(best);
  accuracyResidual(g)=observedAccuracy(g)-predictedAccuracy(g);
  p1=probabilityAtSigma(noiseSDPC(g),B,A,Rpre,groupAdaptation,cfg,z1,z2);
  pRDK1(idx)=p1;
  pChosen=choice1.*p1+(1-choice1).*(1-p1);
  logLikelihood(g)=sum(log(pChosen));
  atBoundary=noiseSDPC(g)<=cfg.noiseBoundsPC(1)*1.002 || ...
    noiseSDPC(g)>=cfg.noiseBoundsPC(2)/1.002;
  if abs(accuracyResidual(g))<=p.Results.MatchTolerance
    matchStatus(g)="matched";
  elseif atBoundary
    matchStatus(g)="boundary";
  else
    matchStatus(g)="closest";
  end

end

groupFit=groupKeys;
groupFit.nTrials=nTrials;
groupFit.observedAccuracy=observedAccuracy;
groupFit.predictedAccuracy=predictedAccuracy;
groupFit.accuracyResidual=accuracyResidual;
groupFit.noiseSDPC=noiseSDPC;
groupFit.logLikelihood=logLikelihood;
groupFit.minGridAccuracy=minGridAccuracy;
groupFit.maxGridAccuracy=maxGridAccuracy;
groupFit.matchStatus=matchStatus;
groupFit.evalCount=evalCount;

fit=struct();
fit.modelVersion=cfg.modelVersion+"-performance-matched";
fit.baseModelVersion=cfg.modelVersion;
selectedAdaptation=unique(adaptationByTrial(selectedIdx));
if numel(selectedAdaptation)==1, fit.adaptation=selectedAdaptation;
else, fit.adaptation=NaN; end
fit.adaptationExponent=cfg.adaptationExponent;
fit.adaptationGainStrength=cfg.adaptationGainStrength;
fit.adaptationReferencePC=cfg.adaptationReferencePC;
fit.groupVariables=string(groupVars);
fit.groupFit=groupFit;
fit.logLikelihood=sum(logLikelihood);
fit.nTrials=numel(selectedIdx);
fit.trialMask=mask;
fit.nNoiseDraws=p.Results.NDraws;
fit.matchTolerance=p.Results.MatchTolerance;
fit.config=cfg;
fit.trialAdaptation=adaptationByTrial;
fit.createdAt=datetime('now');

prediction=struct('pRDK1',pRDK1,'trialIdx',selectedIdx, ...
  'adaptation',fit.adaptation,'groupVariables',string(groupVars), ...
  'adaptationExponent',cfg.adaptationExponent, ...
  'adaptationGainStrength',cfg.adaptationGainStrength, ...
  'adaptationReferencePC',cfg.adaptationReferencePC, ...
  'nNoiseDraws',p.Results.NDraws,'trialAdaptation',adaptationByTrial, ...
  'createdAt',datetime('now'));
end

function acc=accuracyAtSigma(sigma,B,A,Rpre,change1,adaptation,cfg,z1,z2)
pp=probabilityAtSigma(sigma,B,A,Rpre,adaptation,cfg,z1,z2);
acc=mean(change1.*pp+(1-change1).*(1-pp));
end

function pp=probabilityAtSigma(sigma,B,A,Rpre,adaptation,cfg,z1,z2)
pp=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,Rpre,adaptation,sigma,cfg,z1,z2,false);
end
