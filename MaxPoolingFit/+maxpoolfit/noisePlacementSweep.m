function diagnostics = noisePlacementSweep(cache,fit,varargin)
%NOISEPLACEMENTSWEEP Separate pre-max noise from post-selection flips.
%
% For one exact domain-animal-offset-step condition, sweep the Gaussian
% mechanism noise applied before hard-max selection. At every upstream noise
% value, solve for the symmetric categorical flip rate required to reproduce
% observed accuracy, then refit changed/unchanged drift and probe coefficients
% from both the upstream and post-flip choice probabilities.

p=inputParser;
p.addParameter('TrialMask',fit.trialMask,@(x)islogical(x)||isnumeric(x));
p.addParameter('Domain',"IDR",@(x)isstring(x)||ischar(x));
p.addParameter('Animal',"Meetz",@(x)isstring(x)||ischar(x));
p.addParameter('ProbeOffsetDeg',25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('StepMagnitudePC',23,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('NoiseSDGridPC',0:1:24, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('RDKNoiseCorrelation',0, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
p.addParameter('UsePairedRDKDraws',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('MinimumTrials',500, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('MinimumSessions',5, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.parse(varargin{:});

maxpoolfit.validateCache(cache);
assert(isfield(fit,'config')&&isfield(fit,'adaptation')&&isfield(fit,'trialMask'), ...
  'A completed primary fit with configuration and trial mask is required.');
cfg=fit.config; cfg.verbose=false;
T=cache.trials; baseMask=makeMask(T,p.Results.TrialMask);
mask=baseMask & string(T.domain)==string(p.Results.Domain) & ...
  string(T.animal)==string(p.Results.Animal) & ...
  T.probeOffsetDeg==double(p.Results.ProbeOffsetDeg) & ...
  T.stepMagnitudePC==double(p.Results.StepMagnitudePC);
idx=find(mask); assert(numel(idx)>=p.Results.MinimumTrials, ...
  'Only %d selected trials; at least %d are required.',numel(idx),p.Results.MinimumTrials);
nSessions=numel(unique(string(T.sessionID(idx))));
assert(nSessions>=p.Results.MinimumSessions, ...
  'Only %d selected sessions; at least %d are required.',nSessions,p.Results.MinimumSessions);

currentNoise=currentDomainNoise(fit,string(p.Results.Domain));
pCurrent=choiceProbability(cache,idx,fit.adaptation,currentNoise,cfg,p.Results);
[observedGain,currentGain]=gainRows(cache,mask,idx,pCurrent, ...
  p.Results.MinimumTrials,p.Results.MinimumSessions);
observedAccuracy=observedGain.accuracy;
targets=targetTable(observedGain);

noiseGrid=unique(double(p.Results.NoiseSDGridPC(:)),'stable'); n=numel(noiseGrid);
rows=repmat(emptyRow(),n,1);
for k=1:n
  sigma=noiseGrid(k);
  pUp=choiceProbability(cache,idx,fit.adaptation,sigma,cfg,p.Results);
  [~,upGain]=gainRows(cache,mask,idx,pUp, ...
    p.Results.MinimumTrials,p.Results.MinimumSessions);
  upstreamAccuracy=predictedAccuracy(T(idx,:),pUp);
  upstreamLogLikelihood=choiceLogLikelihood(T(idx,:),pUp);
  [epsilon,feasible]=requiredFlip(upstreamAccuracy,observedAccuracy);
  row=emptyRow(); row.noiseSDPC=sigma;
  row.rdkNoiseCorrelation=p.Results.RDKNoiseCorrelation;
  row.nTrials=numel(idx);
  row.nSessions=nSessions; row.observedAccuracy=observedAccuracy;
  row.upstreamAccuracy=upstreamAccuracy; row.requiredFlipRate=epsilon;
  row.downstreamFeasible=feasible;
  row.upstreamLogLikelihood=upstreamLogLikelihood;
  row.upstreamLogLoss=-upstreamLogLikelihood/numel(idx);
  row=assignGain(row,upGain,"upstream");
  if feasible
    pPost=epsilon+(1-2*epsilon).*pUp;
    pPost=min(max(pPost,cfg.probabilityFloor),1-cfg.probabilityFloor);
    [~,postGain]=gainRows(cache,mask,idx,pPost, ...
      p.Results.MinimumTrials,p.Results.MinimumSessions);
    row.postFlipAccuracy=predictedAccuracy(T(idx,:),pPost);
    row.postFlipLogLikelihood=choiceLogLikelihood(T(idx,:),pPost);
    row.postFlipLogLoss=-row.postFlipLogLikelihood/numel(idx);
    row=assignGain(row,postGain,"postFlip");
  end
  row.upstreamRatioRMSE=ratioRMSE(row,targets,"upstream");
  row.postFlipRatioRMSE=ratioRMSE(row,targets,"postFlip");
  if row.downstreamFeasible
    row.ratioRMSE=row.postFlipRatioRMSE;
  else
    row.ratioRMSE=row.upstreamRatioRMSE;
  end
  row.postFlipAbsoluteGainRMSE=absoluteGainRMSE(row,targets);
  rows(k)=row;
end
sweep=struct2table(rows);

finiteRatio=isfinite(sweep.ratioRMSE);
assert(any(finiteRatio),'No sweep value produced finite drift and probe gain ratios.');
objective=sweep.ratioRMSE; objective(~finiteRatio)=inf;
[~,q]=min(objective);
bestRatioRow=sweep(q,:);
feasibleRatio=finiteRatio&sweep.downstreamFeasible;
if any(feasibleRatio)
  objective=sweep.ratioRMSE; objective(~feasibleRatio)=inf;
  [~,q]=min(objective); bestFeasibleRatioRow=sweep(q,:);
else
  bestFeasibleRatioRow=sweep([],:);
end

selection=struct('domain',string(p.Results.Domain), ...
  'animal',string(p.Results.Animal), ...
  'probeOffsetDeg',double(p.Results.ProbeOffsetDeg), ...
  'stepMagnitudePC',double(p.Results.StepMagnitudePC), ...
  'rdkNoiseCorrelation',double(p.Results.RDKNoiseCorrelation), ...
  'pairedRDKDraws',logical(p.Results.UsePairedRDKDraws)|| ...
    p.Results.RDKNoiseCorrelation>0, ...
  'nTrials',numel(idx),'nSessions',nSessions);
diagnostics=struct('sweep',sweep,'targets',targets,'selection',selection, ...
  'trialIdx',idx,'observedAccuracy',observedAccuracy, ...
  'currentNoiseSDPC',currentNoise,'currentGainRow',currentGain, ...
  'bestRatioRow',bestRatioRow,'bestFeasibleRatioRow',bestFeasibleRatioRow, ...
  'adaptation',fit.adaptation,'nNoiseDraws',p.Results.NDraws, ...
  'noiseSDGridPC',noiseGrid,'modelVersion',fit.modelVersion, ...
  'postSelectionRule',"symmetric categorical flip", ...
  'flipEquation',"pFinal = epsilon + (1 - 2*epsilon)*pUpstream", ...
  'createdAt',datetime('now'));
end

function p1=choiceProbability(cache,idx,adaptation,noiseSD,cfg,options)
if options.UsePairedRDKDraws||options.RDKNoiseCorrelation>0
  p1=maxpoolfit.choiceProbabilityCorrelatedRDK(cache,idx,adaptation, ...
    noiseSD,options.RDKNoiseCorrelation,cfg,'NDraws',options.NDraws);
else
  p1=maxpoolfit.choiceProbability(cache,idx,adaptation,noiseSD,cfg, ...
    'NDraws',options.NDraws);
end
end

function mask=makeMask(T,input)
if islogical(input)
  assert(numel(input)==height(T),'Logical TrialMask must contain one value per cached trial.');
  mask=input(:);
else
  mask=false(height(T),1); mask(input)=true;
end
end

function noise=currentDomainNoise(fit,domain)
if isfield(fit,'domainFit')&&istable(fit.domainFit)
  q=string(fit.domainFit.domain)==domain;
  assert(sum(q)==1,'Fit does not contain exactly one noise estimate for %s.',domain);
  noise=double(fit.domainFit.noiseSDPC(q));
else
  q=string(fit.domains)==domain;
  assert(sum(q)==1,'Fit does not contain exactly one noise estimate for %s.',domain);
  noise=double(fit.noiseSDPC(q));
end
end

function [observed,model]=gainRows(cache,mask,idx,p1,minimumTrials,minimumSessions)
prediction=struct('pRDK1',nan(height(cache.trials),1),'trialIdx',idx);
prediction.pRDK1(idx)=p1;
D=maxpoolfit.noiseCoefficientsByStep(cache,prediction, ...
  'TrialMask',mask,'MinimumTrials',minimumTrials, ...
  'ReliableSessionCount',minimumSessions);
S=D.summary; observed=S(string(S.source)=="observed",:);
model=S(string(S.source)=="primary domain noise",:);
assert(height(observed)==1&&height(model)==1, ...
  'Noise-placement selection must resolve to one exact diagnostic cell.');
end

function targets=targetTable(G)
component=["drift";"probe"];
changedCoefficientPerPC=[G.changedDriftCoefficientPerPC;G.changedProbeCoefficientPerPC];
unchangedCoefficientPerPC=[-G.unchangedDriftCoefficientPerPC; ...
  -G.unchangedProbeCoefficientPerPC];
unchangedOverChanged=[G.driftUnchangedOverChanged;G.probeUnchangedOverChanged];
targets=table(component,changedCoefficientPerPC,unchangedCoefficientPerPC, ...
  unchangedOverChanged);
end

function [epsilon,feasible]=requiredFlip(upstreamAccuracy,observedAccuracy)
tolerance=1e-10; epsilon=nan; feasible=false;
if upstreamAccuracy<observedAccuracy-tolerance||upstreamAccuracy<=0.5+tolerance
  return;
end
if abs(upstreamAccuracy-observedAccuracy)<=tolerance
  epsilon=0; feasible=true; return;
end
epsilon=(upstreamAccuracy-observedAccuracy)/(2*upstreamAccuracy-1);
feasible=isfinite(epsilon)&&epsilon>=0&&epsilon<=0.5;
if ~feasible, epsilon=nan; end
end

function accuracy=predictedAccuracy(T,p1)
change1=double(T.changeRDK==1);
accuracy=mean(change1.*p1+(1-change1).*(1-p1));
end

function value=choiceLogLikelihood(T,p1)
y=double(T.choiceRDK==1); pChosen=y.*p1+(1-y).*(1-p1);
pChosen=min(max(pChosen,1e-12),1-1e-12); value=sum(log(pChosen));
end

function row=assignGain(row,G,prefix)
row.(prefix+"ChangedDriftCoefficientPerPC")=G.changedDriftCoefficientPerPC;
row.(prefix+"UnchangedDriftCoefficientPerPC")=-G.unchangedDriftCoefficientPerPC;
row.(prefix+"DriftUnchangedOverChanged")=G.driftUnchangedOverChanged;
row.(prefix+"ChangedProbeCoefficientPerPC")=G.changedProbeCoefficientPerPC;
row.(prefix+"UnchangedProbeCoefficientPerPC")=-G.unchangedProbeCoefficientPerPC;
row.(prefix+"ProbeUnchangedOverChanged")=G.probeUnchangedOverChanged;
end

function value=ratioRMSE(row,targets,prefix)
predicted=[row.(prefix+"DriftUnchangedOverChanged"); ...
  row.(prefix+"ProbeUnchangedOverChanged")];
target=targets.unchangedOverChanged;
if all(isfinite(predicted))&&all(isfinite(target))
  value=sqrt(mean((predicted-target).^2));
else
  value=nan;
end
end

function value=absoluteGainRMSE(row,targets)
predicted=[row.postFlipChangedDriftCoefficientPerPC; ...
  row.postFlipUnchangedDriftCoefficientPerPC; ...
  row.postFlipChangedProbeCoefficientPerPC; ...
  row.postFlipUnchangedProbeCoefficientPerPC];
target=[targets.changedCoefficientPerPC(1);targets.unchangedCoefficientPerPC(1); ...
  targets.changedCoefficientPerPC(2);targets.unchangedCoefficientPerPC(2)];
if all(isfinite(predicted))&&all(isfinite(target))
  value=sqrt(mean((predicted-target).^2));
else
  value=nan;
end
end

function row=emptyRow()
row=struct('noiseSDPC',nan,'rdkNoiseCorrelation',nan, ...
  'nTrials',nan,'nSessions',nan, ...
  'observedAccuracy',nan,'upstreamAccuracy',nan,'requiredFlipRate',nan, ...
  'downstreamFeasible',false,'postFlipAccuracy',nan, ...
  'upstreamLogLikelihood',nan,'postFlipLogLikelihood',nan, ...
  'upstreamLogLoss',nan,'postFlipLogLoss',nan, ...
  'upstreamChangedDriftCoefficientPerPC',nan, ...
  'upstreamUnchangedDriftCoefficientPerPC',nan, ...
  'upstreamDriftUnchangedOverChanged',nan, ...
  'upstreamChangedProbeCoefficientPerPC',nan, ...
  'upstreamUnchangedProbeCoefficientPerPC',nan, ...
  'upstreamProbeUnchangedOverChanged',nan, ...
  'postFlipChangedDriftCoefficientPerPC',nan, ...
  'postFlipUnchangedDriftCoefficientPerPC',nan, ...
  'postFlipDriftUnchangedOverChanged',nan, ...
  'postFlipChangedProbeCoefficientPerPC',nan, ...
  'postFlipUnchangedProbeCoefficientPerPC',nan, ...
  'postFlipProbeUnchangedOverChanged',nan, ...
  'upstreamRatioRMSE',nan,'postFlipRatioRMSE',nan, ...
  'ratioRMSE',nan,'postFlipAbsoluteGainRMSE',nan);
end
