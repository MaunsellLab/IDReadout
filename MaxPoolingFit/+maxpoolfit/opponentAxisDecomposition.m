function result=opponentAxisDecomposition(cache,jointFit,varargin)
%OPPONENTAXISDECOMPOSITION Trace opponent-probe influence through hard max.
%
% Actual paired-probe trials are compared with probe-noise-neutralized
% trials and with two single-180 controls, each retaining one of the actual
% independent probe traces. Identical internal-noise draws are used for
% every paired full/neutral comparison.

p=inputParser;
p.addParameter('ProbeOffsetDeg',170,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('FitSelections',["gain-optimal" "choice-optimal"], ...
  @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('NDraws',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=2));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
if ischar(p.Results.FitSelections), fitNames=string(p.Results.FitSelections);
else, fitNames=unique(string(p.Results.FitSelections(:)),'stable'); end

idx=find(cache.trials.valid & cache.trials.hasStepNoise & ...
  cache.trials.isIncrement & string(cache.trials.domain)=="IDS" & ...
  abs(double(cache.trials.probeOffsetDeg)-p.Results.ProbeOffsetDeg)<1e-9);
assert(~isempty(idx),'No eligible IDS trials found at probe offset %g deg.', ...
  p.Results.ProbeOffsetDeg);
assert(all(cellfun(@numel,cache.stepCohPC(idx,1))>=3) && ...
  all(cellfun(@numel,cache.stepCohPC(idx,2))>=3), ...
  'Opponent-axis decomposition requires paired-probe IDS trials.');

actualCache=cache;
pairedNeutralCache=neutralizeProbes(cache,idx);
singleA=singleProbeCache(cache,idx,2,false);
singleANeutral=singleProbeCache(cache,idx,2,true);
singleB=singleProbeCache(cache,idx,3,false);
singleBNeutral=singleProbeCache(cache,idx,3,true);

summary=table(); trialSummary=table(); gainSummary=table();
for f=1:numel(fitNames)
  [fit,fitLabel]=selectedFit(jointFit,fitNames(f));
  cfg=fit.config; nDraws=fit.nNoiseDraws;
  if ~isempty(p.Results.NDraws), nDraws=p.Results.NDraws; end
  groupRow=find(string(fit.groupFit.analysisGroup)=="IDS",1);
  assert(~isempty(groupRow),'Selected fit does not contain the IDS group.');
  adaptation=fit.groupFit.adaptation(groupRow);
  noiseSD=fit.groupFit.noiseSDPC(groupRow);

  variants={actualCache,pairedNeutralCache,singleA,singleANeutral, ...
    singleB,singleBNeutral};
  labels=["paired actual" "paired neutral" "single 180 A" ...
    "single 180 A neutral" "single 180 B" "single 180 B neutral"];
  evaluations=cell(numel(variants),1);
  for v=1:numel(variants)
    evaluations{v}=evaluateVariant(variants{v},idx,adaptation,noiseSD,cfg,nDraws);
  end

  pairs=[1 2;3 4;5 6];
  pairGain=nan(3,1);
  for q=1:3
    fullIndex=pairs(q,1); neutralIndex=pairs(q,2);
    E=evaluations{fullIndex}; N=evaluations{neutralIndex};
    switchMetrics=winnerSwitches(E,N,actualCache,idx,cfg);
    rows=variantRows(E,N,switchMetrics,fitLabel,labels(fullIndex), ...
      adaptation,noiseSD,nDraws);
    summary=appendRows(summary,rows.summary);
    trialSummary=appendRows(trialSummary,rows.trials);

    prediction=fullPrediction(variants{fullIndex},fit,cfg,nDraws);
    workCache=variants{fullIndex}; workCache.trials.domain=fit.trialGroup;
    C=maxpoolfit.gainComparison(workCache,prediction,prediction, ...
      'TrialMask',fit.trialMask);
    Q=C.summary;
    G=Q(string(Q.source)=="primary domain noise" & ...
      string(Q.domain)=="IDS" & ...
      abs(double(Q.probeOffsetDeg)- ...
      double(workCache.trials.probeOffsetDeg(idx(1))))<1e-9,:);
    assert(~isempty(G),'Predicted opponent-axis gain row was not produced.');
    G.fitSelection=repmat(fitLabel,height(G),1);
    G.configuration=repmat(labels(fullIndex),height(G),1);
    gainSummary=appendRows(gainSummary,G);
    pairGain(q)=mean(G.changedProbeGain,'omitnan');
    target=summary.fitSelection==fitLabel & ...
      summary.configuration==labels(fullIndex);
    summary.predictedChangedProbeGain(target)=pairGain(q);
  end
end

result=struct('summary',summary,'trialSummary',trialSummary, ...
  'gainSummary',gainSummary,'probeOffsetDeg',p.Results.ProbeOffsetDeg, ...
  'trialIdx',idx,'nTrials',numel(idx),'fitSelections',fitNames, ...
  'singleProbeDefinition', ...
    "one actual probe trace retained and placed exactly 180 deg from drift", ...
  'neutralDefinition', ...
    "step-period probe coherence fixed to its inferred pedestal; prestep retained", ...
  'createdAt',datetime('now'));
end

function E=evaluateVariant(cache,idx,a,sigma,cfg,nDraws)
[B,A,P]=maxpoolfit.deterministicResponse(cache,idx,cfg);
D=maxpoolfit.decisionResponse(B,A,P,a,cfg);
[p1,details]=maxpoolfit.choiceProbability(cache,idx,a,sigma,cfg, ...
  'NDraws',nDraws,'ReturnWinners',true);
n=numel(idx); driftResponse=nan(n,1); opponentResponse=nan(n,1);
deterministicProbeWinner=false(n,1); internalProbeWinnerFraction=nan(n,1);
pCorrect=nan(n,1);
for k=1:n
  t=idx(k); changed=double(cache.trials.changeRDK(t));
  drift=double(cache.driftDirDeg(t,changed)); opponent=mod(drift+180,360);
  [~,di]=min(circularDistance(cfg.directionGridDeg,drift));
  [~,oi]=min(circularDistance(cfg.directionGridDeg,opponent));
  driftResponse(k)=D(di,changed,k); opponentResponse(k)=D(oi,changed,k);
  [~,dw]=max(D(:,changed,k));
  deterministicProbeWinner(k)=isProbeWinner(cache,t,changed,dw,cfg);
  if changed==1, winners=details.winnerMechanismIndexRDK1(k,:);
  else, winners=details.winnerMechanismIndexRDK2(k,:); end
  internalProbeWinnerFraction(k)=mean(arrayfun( ...
    @(w)isProbeWinner(cache,t,changed,w,cfg),double(winners)));
  if changed==1, pCorrect(k)=p1(k); else, pCorrect(k)=1-p1(k); end
end
E=struct('cache',cache,'details',details,'pCorrect',pCorrect, ...
  'driftResponse',driftResponse,'opponentResponse',opponentResponse, ...
  'deterministicProbeWinner',deterministicProbeWinner, ...
  'internalProbeWinnerFraction',internalProbeWinnerFraction);
end

function S=winnerSwitches(E,N,cache,idx,cfg)
n=numel(idx); changedMechanism=nan(n,1); driftToProbe=nan(n,1);
probeToDrift=nan(n,1);
for k=1:n
  t=idx(k); changed=double(cache.trials.changeRDK(t));
  if changed==1
    wf=double(E.details.winnerMechanismIndexRDK1(k,:));
    wn=double(N.details.winnerMechanismIndexRDK1(k,:));
  else
    wf=double(E.details.winnerMechanismIndexRDK2(k,:));
    wn=double(N.details.winnerMechanismIndexRDK2(k,:));
  end
  pf=arrayfun(@(w)isProbeWinner(E.cache,t,changed,w,cfg),wf);
  pn=arrayfun(@(w)isProbeWinner(N.cache,t,changed,w,cfg),wn);
  changedMechanism(k)=mean(wf~=wn);
  driftToProbe(k)=mean(~pn&pf); probeToDrift(k)=mean(pn&~pf);
end
S=struct('changedMechanism',changedMechanism, ...
  'driftToProbe',driftToProbe,'probeToDrift',probeToDrift);
end

function R=variantRows(E,N,S,fitLabel,configuration,a,sigma,nDraws)
n=numel(E.pCorrect); probeNoiseDelta=E.pCorrect-N.pCorrect;
summary=table(fitLabel,configuration,n,a,sigma,nDraws, ...
  mean(E.driftResponse),mean(E.opponentResponse), ...
  mean(E.driftResponse-E.opponentResponse), ...
  mean(E.deterministicProbeWinner),mean(E.internalProbeWinnerFraction), ...
  mean(N.internalProbeWinnerFraction),mean(S.changedMechanism), ...
  mean(S.driftToProbe),mean(S.probeToDrift),mean(E.pCorrect), ...
  mean(N.pCorrect),mean(probeNoiseDelta),mean(abs(probeNoiseDelta)),NaN, ...
  'VariableNames',{'fitSelection','configuration','nTrials','adaptation', ...
  'noiseSDPC','nNoiseDraws','meanDriftMechanismResponse', ...
  'meanOpponentMechanismResponse','meanDriftMinusOpponentResponse', ...
  'deterministicProbeWinnerFraction','internalProbeWinnerFraction', ...
  'neutralInternalProbeWinnerFraction','probeNoiseChangesMechanismFraction', ...
  'driftToProbeSwitchFraction','probeToDriftSwitchFraction', ...
  'predictedAccuracy','neutralPredictedAccuracy','meanProbeNoiseDeltaPCorrect', ...
  'meanAbsoluteProbeNoiseDeltaPCorrect','predictedChangedProbeGain'});
trials=table(repmat(fitLabel,n,1),repmat(configuration,n,1), ...
  E.driftResponse,E.opponentResponse,E.driftResponse-E.opponentResponse, ...
  E.internalProbeWinnerFraction,N.internalProbeWinnerFraction, ...
  S.changedMechanism,S.driftToProbe,S.probeToDrift,E.pCorrect,N.pCorrect, ...
  probeNoiseDelta,'VariableNames',{'fitSelection','configuration', ...
  'driftMechanismResponse','opponentMechanismResponse', ...
  'driftMinusOpponentResponse','internalProbeWinnerFraction', ...
  'neutralInternalProbeWinnerFraction','probeNoiseChangesMechanismFraction', ...
  'driftToProbeSwitchFraction','probeToDriftSwitchFraction', ...
  'predictedPCorrect','neutralPredictedPCorrect','probeNoiseDeltaPCorrect'});
R=struct('summary',summary,'trials',trials);
end

function prediction=fullPrediction(cache,fit,cfg,nDraws)
prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg, ...
  'TrialMask',fit.trialMask,'NDraws',nDraws);
end

function cache=neutralizeProbes(cache,idx)
for q=1:numel(idx)
  t=idx(q);
  for rdk=1:2
    values=double(cache.stepCohPC{t,rdk}(:));
    for stream=2:numel(values)
      values(stream)=probePedestal(cache,t,rdk,stream);
    end
    cache.stepCohPC{t,rdk}=values;
  end
end
end

function cache=singleProbeCache(cache,idx,stream,neutral)
for q=1:numel(idx)
  t=idx(q);
  for rdk=1:2
    drift=double(cache.driftDirDeg(t,rdk));
    step=double(cache.stepCohPC{t,rdk}(:));
    pre=double(cache.preCohPC{t,rdk}(:));
    probeValue=step(stream);
    if neutral, probeValue=probePedestal(cache,t,rdk,stream); end
    cache.streamDirDeg{t,rdk}=[drift;mod(drift+180,360)];
    cache.stepCohPC{t,rdk}=[step(1);probeValue];
    cache.preCohPC{t,rdk}=[pre(1);pre(stream)];
  end
  cache.trials.probeOffsetDeg(t)=180;
end
end

function pedestal=probePedestal(cache,t,rdk,stream)
pre=double(cache.preCohPC{t,rdk}(:));
pedestal=pre(stream);
if isfield(cache,'raw')&&isfield(cache.raw,'stepStreamCohPC')
  X=double(cache.raw.stepStreamCohPC{t,rdk}(:,stream));
  X=X(isfinite(X));
  if ~isempty(X)&&max(X)>min(X)
    pedestal=.5*(min(X)+max(X));
  end
end
end

function tf=isProbeWinner(cache,t,rdk,winnerIndex,cfg)
winnerDir=double(cfg.directionGridDeg(double(winnerIndex)));
drift=double(cache.driftDirDeg(t,rdk));
streams=double(cache.streamDirDeg{t,rdk}(:));
probe=streams(circularDistance(streams,drift)>1e-6);
if isempty(probe), tf=false; return; end
probeDistance=min(circularDistance(winnerDir,probe));
tf=probeDistance<circularDistance(winnerDir,drift);
end

function d=circularDistance(a,b)
d=abs(mod(double(a)-double(b)+180,360)-180);
end

function [fit,label]=selectedFit(result,name)
switch string(name)
  case "gain-optimal", fit=result.gainOptimalFit; label="gain-optimal";
  case "choice-optimal", fit=result.choiceOptimalFit; label="choice-optimal";
  otherwise, error('Unknown FitSelections value %s.',string(name));
end
assert(isfield(fit,'adaptationByGroup')&&isfield(fit,'trialGroup'), ...
  'Selected joint-fit result is incomplete.');
end

function A=appendRows(A,B)
if width(A)==0, A=B; else, A=[A;B]; end
end
