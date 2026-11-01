function diagnostics = noiseCoefficientsByStep(cache,prediction,varargin)
%NOISECOEFFICIENTSBYSTEP Estimate RDK-specific noise effects at exact steps.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('MinimumTrials',100, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('ReliableSessionCount',5, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
T=cache.trials;
if isempty(p.Results.TrialMask)
  mask=T.valid & T.hasStepNoise & T.isIncrement;
elseif islogical(p.Results.TrialMask)
  assert(numel(p.Results.TrialMask)==height(T), ...
    'Logical TrialMask must contain one value per cached trial.');
  mask=p.Results.TrialMask(:);
else
  mask=false(height(T),1); mask(p.Results.TrialMask)=true;
end
idx=find(mask);
assert(all(isfinite(prediction.pRDK1(idx))), ...
  'Prediction does not cover every selected trial.');
P=maxpoolfit.gainPredictors(cache,idx);
observedCorrect=double(T.choiceRDK(idx)==T.changeRDK(idx));
change1=double(T.changeRDK(idx)==1);
predictedCorrect=change1.*double(prediction.pRDK1(idx))+ ...
  (1-change1).*(1-double(prediction.pRDK1(idx)));
outcomes=[observedCorrect predictedCorrect];
sources=["observed" "primary domain noise"];

groupVars={'domain','animal','probeOffsetDeg','stepMagnitudePC'};
groups=unique(P(:,groupVars),'rows','stable');
summary=emptySummary(); fits=cell(height(groups),numel(sources));
for g=1:height(groups)
  q=string(P.domain)==string(groups.domain(g)) & ...
    string(P.animal)==string(groups.animal(g)) & ...
    P.probeOffsetDeg==groups.probeOffsetDeg(g) & ...
    P.stepMagnitudePC==groups.stepMagnitudePC(g);
  n=sum(q);
  if n<p.Results.MinimumTrials, continue; end
  Pg=P(q,:); Y=outcomes(q,:);
  [X,coefficientColumns,sessionGroup,nSessions]=designMatrix(Pg);
  for s=1:numel(sources)
    [b,V,fitInfo]=fitBinomial(X,Y(:,s),sessionGroup);
    fits{g,s}=fitInfo;
    row=groups(g,:); row.source=sources(s);
    row.nTrials=n; row.nSessions=nSessions;
    row.inferenceReliable=nSessions>=p.Results.ReliableSessionCount;
    row.accuracy=mean(Y(:,s));
    names=["changedDrift" "unchangedDrift" ...
      "changedProbe" "unchangedProbe"];
    for k=1:numel(names)
      column=coefficientColumns(k); prefix=names(k);
      row.(prefix+"CoefficientPerPC")=b(column);
      row.(prefix+"CoefficientSE")=sqrt(max(V(column,column),0));
    end
    row.driftUnchangedOverChanged= ...
      -row.unchangedDriftCoefficientPerPC/row.changedDriftCoefficientPerPC;
    row.probeUnchangedOverChanged= ...
      -row.unchangedProbeCoefficientPerPC/row.changedProbeCoefficientPerPC;
    [row.driftGainImpliedWinFraction,row.driftGainImpliedWinFractionSE, ...
      row.driftGainSignsConsistent]=gainWinFraction( ...
      b,V,coefficientColumns(1),coefficientColumns(2));
    [row.probeGainImpliedWinFraction,row.probeGainImpliedWinFractionSE, ...
      row.probeGainSignsConsistent]=gainWinFraction( ...
      b,V,coefficientColumns(3),coefficientColumns(4));
    if isempty(summary), summary=row; else, summary=[summary;row]; end %#ok<AGROW>
  end
end

diagnostics=struct('summary',summary,'sources',sources, ...
  'trialIdx',idx,'nTrials',numel(idx), ...
  'minimumTrials',p.Results.MinimumTrials, ...
  'reliableSessionCount',p.Results.ReliableSessionCount, ...
  'coefficientDefinition', ...
    "binomial coefficient per percentage-point kernel-weighted coherence", ...
  'regressionDefinition', ...
    "exact step and offset; session fixed effects; within-session centered noise predictors", ...
  'fitInfo',{fits},'createdAt',datetime('now'));
end

function T=emptySummary()
names={'domain','animal','probeOffsetDeg','stepMagnitudePC','source', ...
  'nTrials','nSessions','inferenceReliable','accuracy', ...
  'changedDriftCoefficientPerPC','changedDriftCoefficientSE', ...
  'unchangedDriftCoefficientPerPC','unchangedDriftCoefficientSE', ...
  'changedProbeCoefficientPerPC','changedProbeCoefficientSE', ...
  'unchangedProbeCoefficientPerPC','unchangedProbeCoefficientSE', ...
  'driftUnchangedOverChanged','probeUnchangedOverChanged', ...
  'driftGainImpliedWinFraction','driftGainImpliedWinFractionSE', ...
  'driftGainSignsConsistent', ...
  'probeGainImpliedWinFraction','probeGainImpliedWinFractionSE', ...
  'probeGainSignsConsistent'};
types=[{'string','string','double','double','string','double','double', ...
  'logical'} repmat({'double'},1,13) {'logical'} ...
  repmat({'double'},1,2) {'logical'}];
T=table('Size',[0 numel(names)],'VariableTypes',types, ...
  'VariableNames',names);
end

function [w,se,consistent]=gainWinFraction(b,V,changedColumn,unchangedColumn)
c=b(changedColumn); u=b(unchangedColumn);
consistent=isfinite(c)&&isfinite(u)&&c>0&&u<0;
w=nan; se=nan;
if ~consistent, return; end
denominator=c-u;
w=c/denominator;
gradient=[-u c]./(denominator.^2);
covariance=V([changedColumn unchangedColumn], ...
  [changedColumn unchangedColumn]);
if all(isfinite(covariance),'all')
  variance=gradient*covariance*gradient';
  se=sqrt(max(variance,0));
end
end

function [X,columns,sessionGroup,nSessions]=designMatrix(P)
n=height(P); [sessionGroup,sessions]=findgroups(string(P.sessionID));
nSessions=numel(sessions);
sessionIntercept=sparse((1:n)',sessionGroup,1,n,nSessions);
sideSign=2.*double(P.changeRDK==1)-1;
noise=[P.changedDriftNoisePC P.unchangedDriftNoisePC ...
  P.changedProbeNoisePC P.unchangedProbeNoisePC];
for k=1:4, noise(:,k)=centerWithinGroup(noise(:,k),sessionGroup); end
X=[sessionIntercept sparse(sideSign) sparse(noise)];
columns=(size(X,2)-3):size(X,2);
end

function centered=centerWithinGroup(x,group)
means=splitapply(@mean,x,group);
centered=x-means(group);
end

function [b,V,info]=fitBinomial(X,y,cluster)
y=double(y(:)); assert(all(isfinite(y)&y>=0&y<=1));
[n,nParameters]=size(X); b=zeros(nParameters,1); ridge=1e-8;
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
H=full(Xw'*Xw)+ridge.*eye(nParameters); bread=pinv(H);
groups=unique(cluster); meat=zeros(nParameters);
for k=1:numel(groups)
  q=cluster==groups(k); score=full(X(q,:)'*(y(q)-mu(q)));
  meat=meat+score*score';
end
if numel(groups)>=2 && n>nParameters
  V=bread*meat*bread;
  V=V.*(numel(groups)/(numel(groups)-1)).*((n-1)/(n-nParameters));
else
  V=nan(nParameters);
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
