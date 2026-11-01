function comparison = gainComparison(cache,primaryPrediction,matchedPrediction,varargin)
%GAINCOMPARISON Estimate observed and predicted choice-noise gains identically.
%
% Noise coefficients are divided by an observed empirical psychometric slope,
% yielding equivalent step-coherence gain. Session fixed effects remove
% pedestal constants and day-to-day threshold differences. Changed- and
% unchanged-RDK predictors remain separate.

p=inputParser;
p.addParameter('TrialMask',[],@(x)isempty(x)||islogical(x)||isnumeric(x));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
T=cache.trials;
required={'stepMagnitudePC','probeOffsetDeg','baseCohPC','stepEndpointCohPC'};
assert(all(ismember(required,T.Properties.VariableNames)), ...
  'Cache lacks fields required for physical-coherence gain estimation.');

if isempty(p.Results.TrialMask)
  mask=T.valid & T.hasStepNoise & T.isIncrement;
elseif islogical(p.Results.TrialMask)
  mask=p.Results.TrialMask(:);
else
  mask=false(height(T),1); mask(p.Results.TrialMask)=true;
end
assert(numel(mask)==height(T),'TrialMask length does not match the cache.');
idx=find(mask);
assert(all(isfinite(primaryPrediction.pRDK1(idx))) && ...
  all(isfinite(matchedPrediction.pRDK1(idx))), ...
  'Both predictions must cover every gain-analysis trial.');

P=maxpoolfit.gainPredictors(cache,idx);
observedCorrect=double(T.choiceRDK(idx)==T.changeRDK(idx));
primaryCorrect=probabilityCorrect(T.changeRDK(idx),primaryPrediction.pRDK1(idx));
matchedCorrect=probabilityCorrect(T.changeRDK(idx),matchedPrediction.pRDK1(idx));
outcomes=[observedCorrect primaryCorrect matchedCorrect];
sources=["observed" "primary domain noise" "performance matched"];

groups=unique(P(:,{'domain','animal'}),'rows','stable');
summary=table(); reference=groups;
reference.nTrials=zeros(height(groups),1);
reference.nSessions=zeros(height(groups),1);
reference.slopePerPC=nan(height(groups),1);
reference.slopeSE=nan(height(groups),1);
reference.method=strings(height(groups),1);
fits=cell(height(groups),numel(sources));
for g=1:height(groups)
  groupAll=string(T.domain)==string(groups.domain(g)) & ...
    string(T.animal)==string(groups.animal(g));
  refMask=groupAll & T.valid & ~T.hasStepNoise & T.isIncrement;
  noisyRefMask=groupAll & T.valid & T.hasStepNoise & T.isIncrement;
  [refSlope,refVariance,refInfo]=referenceSlope( ...
    T(refMask,:),T(noisyRefMask,:));
  assert(isfinite(refSlope) && refSlope>0, ...
    'No positive empirical reference slope for %s - %s.', ...
    groups.domain(g),groups.animal(g));
  reference.nTrials(g)=refInfo.nTrials;
  reference.nSessions(g)=refInfo.nSessions;
  reference.slopePerPC(g)=refSlope;
  reference.slopeSE(g)=sqrt(max(refVariance,0));
  reference.method(g)=refInfo.method;

  q=string(P.domain)==string(groups.domain(g)) & ...
    string(P.animal)==string(groups.animal(g));
  Pg=P(q,:); Y=outcomes(q,:);
  [X,offsets,coefficientColumns,designInfo]=gainDesign(Pg);
  for s=1:numel(sources)
    [b,V,fitInfo]=fitBinomial(X,Y(:,s),designInfo.sessionGroup);
    fits{g,s}=fitInfo;
    for k=1:numel(offsets)
      row=groups(g,:);
      row.source=sources(s);
      row.probeOffsetDeg=offsets(k);
      oq=Pg.probeOffsetDeg==offsets(k);
      row.nTrials=sum(oq);
      row.nSessions=numel(unique(Pg.sessionID(oq)));
      if row.nSessions<3
        row.uncertaintyStatus="preliminary: fewer than 3 sessions";
      else
        row.uncertaintyStatus="session-clustered CI available";
      end
      row.referenceSlopePerPC=refSlope;
      names={"changedDrift","unchangedDrift","changedProbe","unchangedProbe"};
      for c=1:4
        col=coefficientColumns(k,c);
        coefficient=b(col);
        coefficientVariance=max(V(col,col),0);
        gain=coefficient/refSlope;
        gainVariance=coefficientVariance/(refSlope^2)+ ...
          (coefficient^2)*refVariance/(refSlope^4);
        prefix=names{c};
        row.(prefix+"CoefficientPerPC")=coefficient;
        row.(prefix+"Gain")=gain;
        row.(prefix+"GainSE")=sqrt(max(gainVariance,0));
        if row.nSessions<2
          % A cluster-robust variance is not identifiable from one session.
          row.(prefix+"GainSE")=nan;
        end
        row.(prefix+"GainCILow")=gain-1.96*row.(prefix+"GainSE");
        row.(prefix+"GainCIHigh")=gain+1.96*row.(prefix+"GainSE");
      end
      row.changedProbeOverDrift=row.changedProbeGain/row.changedDriftGain;
      if isempty(summary), summary=row; else, summary=[summary;row]; end %#ok<AGROW>
    end
  end
end

comparison=struct();
comparison.summary=summary;
comparison.referencePsychometric=reference;
comparison.sources=sources;
comparison.trialIdx=idx;
comparison.nTrials=numel(idx);
comparison.probePredictorContract=unique(P(:,{'domain','probePredictorContract'}),'rows');
comparison.fitInfo=fits;
comparison.minimumSessionsForDisplayedCI=3;
comparison.gainDefinition= ...
  "binomial choice-noise coefficient divided by an observed empirical "+ ...
  "psychometric slope; gains are equivalent step-coherence units";
comparison.regressionDefinition= ...
  "session fixed effects; within-session centered predictors; separate "+ ...
  "changed/unchanged drift and probe coefficients by offset";
comparison.createdAt=datetime('now');
end

function pCorrect=probabilityCorrect(changeRDK,pRDK1)
change1=double(changeRDK==1);
pCorrect=change1.*double(pRDK1)+(1-change1).*(1-double(pRDK1));
end

function [slope,slopeVariance,info]=referenceSlope(noNoiseT,noisyT)
% Prefer the planned no-noise psychometric reference. ProbeSession caches may
% omit those trials, in which case use noisy trials with offset fixed effects.
useNoNoise=height(noNoiseT)>=20;
if useNoNoise
  [sessionGroup,sessions]=findgroups(string(noNoiseT.sessionID));
  step=centerWithinGroup(double(noNoiseT.stepMagnitudePC),sessionGroup);
  useNoNoise=std(step)>1e-6;
end
if useNoNoise
  T=noNoiseT; n=height(T); nSessions=numel(sessions);
  intercept=sparse((1:n)',sessionGroup,1,n,nSessions);
  method="no-step-noise trials with session fixed effects";
else
  T=noisyT;
  assert(height(T)>=20, ...
    'Too few increment trials for an empirical reference slope.');
  [sessionGroup,sessions]=findgroups(string(T.sessionID));
  [offsetGroup,offsets]=findgroups(double(T.probeOffsetDeg));
  step=centerWithinGroup(double(T.stepMagnitudePC),offsetGroup);
  assert(std(step)>1e-6, ...
    'No within-offset step-magnitude variation for an empirical reference slope.');
  n=height(T); nSessions=numel(sessions);
  intercept=sparse((1:n)',offsetGroup,1,n,numel(offsets));
  method="step-noise trials with probe-offset fixed effects";
end
sideSign=2.*double(T.changeRDK==1)-1;
X=[intercept sparse(sideSign) sparse(step)];
y=double(T.choiceRDK==T.changeRDK);
[b,V,fitInfo]=fitBinomial(X,y,sessionGroup);
slope=b(end); slopeVariance=V(end,end);
info=fitInfo; info.nTrials=n; info.nSessions=nSessions; info.method=method;
end

function [X,offsets,coefficientColumns,info]=gainDesign(P)
n=height(P); [sessionGroup,sessions]=findgroups(string(P.sessionID));
offsets=unique(double(P.probeOffsetDeg),'sorted'); nOffsets=numel(offsets);
sessionIntercept=sparse((1:n)',sessionGroup,1,n,numel(sessions));
sideSign=2.*double(P.changeRDK==1)-1;
step=centerWithinGroup(double(P.stepMagnitudePC),sessionGroup);
noise=[P.changedDriftNoisePC P.unchangedDriftNoisePC ...
  P.changedProbeNoisePC P.unchangedProbeNoisePC];
for c=1:size(noise,2), noise(:,c)=centerWithinGroup(noise(:,c),sessionGroup); end
X=[sessionIntercept sparse(sideSign) sparse(step)];
coefficientColumns=zeros(nOffsets,4);
for k=1:nOffsets
  q=double(P.probeOffsetDeg==offsets(k));
  for c=1:4
    X=[X sparse(q.*noise(:,c))]; %#ok<AGROW>
    coefficientColumns(k,c)=size(X,2);
  end
end
info=struct('sessionGroup',sessionGroup,'sessions',sessions, ...
  'nSessions',numel(sessions),'nOffsets',nOffsets);
end

function centered=centerWithinGroup(x,group)
means=splitapply(@mean,x,group);
centered=x-means(group);
end

function [b,V,info]=fitBinomial(X,y,cluster)
% Newton scoring for binary or fractional binomial outcomes.
y=double(y(:)); assert(all(isfinite(y) & y>=0 & y<=1));
[n,nParameters]=size(X); b=zeros(nParameters,1);
ridge=1e-8;
converged=false; iteration=0;
for iteration=1:100
  eta=X*b; mu=logistic(eta); w=max(mu.*(1-mu),1e-6);
  Xw=spdiags(sqrt(w),0,n,n)*X;
  H=full(Xw'*Xw)+ridge.*eye(nParameters);
  delta=H\full(X'*(y-mu));
  oldNLL=binomialNLL(y,mu); scale=1;
  while scale>=1/1024
    candidate=b+scale.*delta;
    if binomialNLL(y,logistic(X*candidate))<=oldNLL, break; end
    scale=scale/2;
  end
  b=candidate;
  if max(abs(scale.*delta))<1e-7, converged=true; break; end
end
eta=X*b; mu=logistic(eta); w=max(mu.*(1-mu),1e-6);
Xw=spdiags(sqrt(w),0,n,n)*X;
H=full(Xw'*Xw)+ridge.*eye(nParameters);
bread=pinv(H); groups=unique(cluster); meat=zeros(nParameters);
for k=1:numel(groups)
  q=cluster==groups(k); score=full(X(q,:)'*(y(q)-mu(q)));
  meat=meat+score*score';
end
V=bread*meat*bread;
if numel(groups)>1 && n>nParameters
  V=V.*(numel(groups)/(numel(groups)-1)).*((n-1)/(n-nParameters));
end
info=struct('converged',converged,'iterations',iteration, ...
  'nTrials',n,'nParameters',nParameters,'nClusters',numel(groups), ...
  'logLikelihood',-binomialNLL(y,mu));
end

function p=logistic(x)
p=1./(1+exp(-max(min(x,35),-35)));
end

function value=binomialNLL(y,p)
p=min(max(p,1e-12),1-1e-12);
value=-sum(y.*log(p)+(1-y).*log(1-p));
end
