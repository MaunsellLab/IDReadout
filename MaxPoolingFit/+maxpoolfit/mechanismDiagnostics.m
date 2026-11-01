function diagnostics = mechanismDiagnostics( ...
    cache,primaryFit,performanceFit,conditionRows,varargin)
%MECHANISMDIAGNOSTICS Audit winners and competitive gaps in selected conditions.

% Direction mechanisms are assigned to the nearest displayed stream center.
% The drift region contains mechanisms nearer the drift stream than any probe
% stream; all remaining mechanisms form the probe region.

p=inputParser;
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MaxTrialsPerCondition',1000, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.addParameter('NoiseGridPC',[0 1 2 3 4 5 6 8 10 12 15 20 30 40], ...
  @(x)isnumeric(x)&&isvector(x)&&all(x>=0));
p.addParameter('Seed',4817,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('WinnerBinEdgesDeg',-180:10:180, ...
  @(x)isnumeric(x)&&isvector(x)&&numel(x)>=3);
p.parse(varargin{:});
cfg=performanceFit.config;
assert(isfinite(cfg.adaptation) && isfinite(primaryFit.adaptation) && ...
  abs(cfg.adaptation-primaryFit.adaptation)<1e-12, ...
  'Mechanism diagnostics require matching shared adaptation values.');
G=performanceFit.groupFit;
assert(all(conditionRows>=1 & conditionRows<=height(G)), ...
  'ConditionRows contains invalid table rows.');
conditionRows=unique(conditionRows(:),'stable');
groupVars=cellstr(performanceFit.groupVariables);
z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
theta=double(cfg.directionGridDeg(:));
edges=double(p.Results.WinnerBinEdgesDeg(:)');
centers=(edges(1:end-1)+edges(2:end))/2;

old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(p.Results.Seed,'twister');
nConditions=numel(conditionRows);
S=G(conditionRows,groupVars);
S.nTrialsTotal=zeros(nConditions,1);
S.nTrialsDiagnostic=zeros(nConditions,1);
S.observedAccuracy=zeros(nConditions,1);
S.primaryPredictedAccuracy=zeros(nConditions,1);
S.performancePredictedAccuracy=zeros(nConditions,1);
S.primaryNoiseSDPC=zeros(nConditions,1);
S.performanceNoiseSDPC=zeros(nConditions,1);
S.diagnosticPredictedAccuracy=zeros(nConditions,1);
S.deterministicAccuracy=zeros(nConditions,1);
S.changedDriftWinFraction=zeros(nConditions,1);
S.unchangedDriftWinFraction=zeros(nConditions,1);
S.bothRDKsProbeWinFraction=zeros(nConditions,1);
S.meanChangedDriftProbeGap=zeros(nConditions,1);
S.meanUnchangedDriftProbeGap=zeros(nConditions,1);
S.meanBetweenRDKMargin=zeros(nConditions,1);
S.medianBetweenRDKMargin=zeros(nConditions,1);

histogram=table(); sweep=table();
for c=1:nConditions
  sourceRow=conditionRows(c);
  idx=find(conditionMask(cache.trials,performanceFit.trialMask, ...
    G(sourceRow,:),groupVars));
  S.nTrialsTotal(c)=numel(idx);
  if numel(idx)>p.Results.MaxTrialsPerCondition
    idx=idx(randperm(numel(idx),p.Results.MaxTrialsPerCondition));
    idx=sort(idx);
  end
  S.nTrialsDiagnostic(c)=numel(idx);
  S.observedAccuracy(c)=G.observedAccuracy(sourceRow);
  S.primaryPredictedAccuracy(c)=G.primaryPredictedAccuracy(sourceRow);
  S.performancePredictedAccuracy(c)=G.predictedAccuracy(sourceRow);
  S.performanceNoiseSDPC(c)=G.noiseSDPC(sourceRow);
  domain=string(cache.trials.domain(idx(1)));
  primaryRow=find(string(primaryFit.domains)==domain,1);
  assert(~isempty(primaryRow),'Primary fit has no noise SD for domain %s.',domain);
  S.primaryNoiseSDPC(c)=primaryFit.noiseSDPC(primaryRow);
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  M=summarizeCondition(cache,idx,B,A,Rpre,S.performanceNoiseSDPC(c), ...
    cfg,z1,z2,theta,edges);
  names=setdiff(M.Properties.VariableNames, ...
    {'winnerBinCenterDeg','changedWinnerFraction','unchangedWinnerFraction'});
  for k=1:numel(names), S.(names{k})(c)=M.(names{k})(1); end

  H=S(c,groupVars);
  nBins=numel(centers);
  H=repmat(H,nBins,1);
  H.winnerBinCenterDeg=centers(:);
  H.changedWinnerFraction=M.changedWinnerFraction(:);
  H.unchangedWinnerFraction=M.unchangedWinnerFraction(:);
  histogram=appendTable(histogram,H);

  sigmas=unique([double(p.Results.NoiseGridPC(:)); ...
    S.primaryNoiseSDPC(c);S.performanceNoiseSDPC(c)]);
  change1=double(cache.trials.changeRDK(idx)==1);
  acc=zeros(numel(sigmas),1);
  for k=1:numel(sigmas)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      B,A,Rpre,cfg.adaptation,sigmas(k),cfg,z1,z2,false);
    acc(k)=mean(change1.*p1+(1-change1).*(1-p1));
  end
  W=S(c,groupVars); W=repmat(W,numel(sigmas),1);
  W.noiseSDPC=sigmas;
  W.predictedAccuracy=acc;
  W.observedAccuracy=repmat(S.observedAccuracy(c),numel(sigmas),1);
  sweep=appendTable(sweep,W);
end

diagnostics=struct('conditionSummary',S,'noiseSweep',sweep, ...
  'winnerHistogram',histogram,'conditionRows',conditionRows, ...
  'nDraws',p.Results.NDraws, ...
  'maxTrialsPerCondition',p.Results.MaxTrialsPerCondition, ...
  'seed',p.Results.Seed,'createdAt',datetime('now'));
end

function M=summarizeCondition(cache,idx,B,A,Rpre,sigma,cfg,z1,z2,theta,edges)
n=numel(idx); nDraws=size(z1,2);
gapChanged=zeros(n,nDraws); gapUnchanged=zeros(n,nDraws);
margin=zeros(n,nDraws); changedDrift=zeros(n,nDraws); unchangedDrift=zeros(n,nDraws);
hChanged=zeros(1,numel(edges)-1); hUnchanged=hChanged;
detCorrect=zeros(n,1);
for q=1:n
  D=maxpoolfit.decisionResponse(B(:,:,q),A(:,:,q),Rpre(:,:,q), ...
    cfg.adaptation,cfg);
  R=cell(2,1); gap=cell(2,1); winner=cell(2,1);
  for j=1:2
    Z=z1; if j==2, Z=z2; end
    R{j}=D(:,j)+sigma.*Z;
    [maxAll,winner{j}]=max(R{j},[],1);
    drift=double(cache.driftDirDeg(idx(q),j));
    streams=double(cache.streamDirDeg{idx(q),j}(:));
    probe=streams(circularDistance(streams,drift)>1e-6);
    assert(~isempty(probe),'Selected diagnostic trial has no probe stream.');
    driftMask=circularDistance(theta,drift)<=min( ...
      circularDistance(theta,probe'),[],2);
    maxDrift=max(R{j}(driftMask,:),[],1);
    maxProbe=max(R{j}(~driftMask,:),[],1);
    gap{j}=maxDrift-maxProbe;
    R{j}=maxAll;
  end
  if cache.trials.changeRDK(idx(q))==1
    ch=1; un=2;
  else
    ch=2; un=1;
  end
  gapChanged(q,:)=gap{ch}; gapUnchanged(q,:)=gap{un};
  changedDrift(q,:)=gap{ch}>=0; unchangedDrift(q,:)=gap{un}>=0;
  margin(q,:)=R{ch}-R{un};
  detMax=[max(D(:,1)) max(D(:,2))];
  detCorrect(q)=double(detMax(ch)>detMax(un))+ ...
    0.5*double(detMax(ch)==detMax(un));
  relChanged=signedCircularDifference(theta(winner{ch}), ...
    cache.driftDirDeg(idx(q),ch));
  relUnchanged=signedCircularDifference(theta(winner{un}), ...
    cache.driftDirDeg(idx(q),un));
  hChanged=hChanged+histcounts(relChanged,edges);
  hUnchanged=hUnchanged+histcounts(relUnchanged,edges);
end
M=table();
M.diagnosticPredictedAccuracy=mean(margin(:)>0)+0.5*mean(margin(:)==0);
M.deterministicAccuracy=mean(detCorrect);
M.changedDriftWinFraction=mean(changedDrift,'all');
M.unchangedDriftWinFraction=mean(unchangedDrift,'all');
M.bothRDKsProbeWinFraction=mean(~changedDrift & ~unchangedDrift,'all');
M.meanChangedDriftProbeGap=mean(gapChanged,'all');
M.meanUnchangedDriftProbeGap=mean(gapUnchanged,'all');
M.meanBetweenRDKMargin=mean(margin,'all');
M.medianBetweenRDKMargin=median(margin,'all');
M.winnerBinCenterDeg=(edges(1:end-1)+edges(2:end))/2;
M.changedWinnerFraction=hChanged/sum(hChanged);
M.unchangedWinnerFraction=hUnchanged/sum(hUnchanged);
end

function mask=conditionMask(T,baseMask,row,groupVars)
mask=baseMask;
for k=1:numel(groupVars)
  name=groupVars{k}; value=row.(name)(1);
  if isnumeric(T.(name))
    tolerance=1e-5*max(1,abs(double(value)));
    mask=mask & abs(double(T.(name))-double(value))<=tolerance;
  else
    mask=mask & string(T.(name))==string(value);
  end
end
end

function d=circularDistance(a,b)
d=abs(mod(double(a)-double(b)+180,360)-180);
end

function d=signedCircularDifference(a,b)
d=mod(double(a)-double(b)+180,360)-180;
end

function A=appendTable(A,B)
if isempty(A), A=B; else, A=[A;B]; end
end
