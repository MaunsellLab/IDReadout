function F = analyzeIDRControlledPooling( ...
  D, signedGainFits, betaWeibull, lapse, targetPerformance, gainBounds, varargin)
% analyzeIDRControlledPooling  Compare IDR sum, p-norm, and max pooling.
%
% This is a controlled empirical comparison. The change-side contribution
% is represented flexibly by four signed gains in every model. Only the
% no-change pooling rule differs.
%
% For each no-change stream x, opponent-consistent candidate magnitude is
%
%   z(x) = max(x, -x/r),
%
% with fixed preferred:null ratio r. Preferred and probe candidate scales
% are fit separately and constrained nonnegative. For paired probe offsets,
% q is the single-candidate predictor and nProbeCandidates=2:
%
%   sum:    aP*zP + nProbeCandidates*aQ*zQ
%   p-norm: ((aP*zP)^p + nProbeCandidates*(aQ*zQ)^p)^(1/p)
%   max:    max(aP*zP, aQ*zQ)
%
% Name-value arguments:
%   PreferredNullRatio  Default 2.42
%   PBounds             Bounds for fitted p (default [1 100])
%   ProfilePValues      Fixed-p likelihood-profile grid
%   PlotPath            Optional summary PDF path


p=inputParser;
addParameter(p,'PreferredNullRatio',2.42, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0);
addParameter(p,'PBounds',[1 100], ...
  @(x) isnumeric(x) && numel(x)==2 && all(isfinite(x)) && x(1)>=1 && x(2)>x(1));
addParameter(p,'ProfilePValues',[1 1.5 2 3 5 10 20 50 100], ...
  @(x) isnumeric(x) && isvector(x) && ~isempty(x) && all(isfinite(x)) && all(x>=1));
addParameter(p,'PlotPath','',@(x) ischar(x) || isstring(x));
parse(p,varargin{:});
opts=p.Results;
profilePValues=unique(double(opts.ProfilePValues(:)),'sorted');
if any(profilePValues<opts.PBounds(1) | profilePValues>opts.PBounds(2))
  error('analyzeIDRControlledPooling:ProfileOutsideBounds', ...
    'Every ProfilePValues entry must lie within PBounds.');
end

validateInput(D,signedGainFits,betaWeibull,lapse,targetPerformance,gainBounds);

T=D.trialTable;
stepFrames=D.stepFrames(:);
CP=stepMean(D.changePrefNoiseByFrameTrial,stepFrames);
CQ=stepMean(D.changeProbeEffectiveNoiseByFrameTrial,stepFrames);
NP=stepMean(D.noChangePrefNoiseByFrameTrial,stepFrames);
NQcandidate=stepMean(D.noChangeProbeCandidateNoiseByFrameTrial,stepFrames);

XchangeAll=[posPart(CP),negPart(CP),posPart(CQ),negPart(CQ)];
r=opts.PreferredNullRatio;
zPAll=max(NP,-NP/r);
zQAll=max(NQcandidate,-NQcandidate/r);

offsets=unique(double(T.probeDirDeg),'sorted');
offsets=offsets(isfinite(offsets));
nOffsets=numel(offsets);
emptyResult=struct('probeDirDeg',[],'nTrials',[],'nSessions',[], ...
  'nProbeCandidates',[],'sum',[],'pNorm',[],'hardMax',[], ...
  'freeSignedNLL',[],'sessionIndex',[]);
byOffset=repmat(emptyResult,nOffsets,1);

nAll=height(T);
trialNLLSum=nan(nAll,1); trialNLLPNorm=nan(nAll,1); trialNLLHardMax=nan(nAll,1);

for k=1:nOffsets
  theta=offsets(k);
  use=double(T.probeDirDeg)==theta;
  nProbeCandidates=unique(double(T.nYokedProbeStreams(use)));
  if numel(nProbeCandidates)~=1
    error('analyzeIDRControlledPooling:MixedProbeCount', ...
      'Offset %g has inconsistent probe-candidate counts.',theta);
  end

  common=struct();
  common.stepCoh=double(T.stepDeltaCohPC(use));
  common.correct=double(T.correct(use));
  common.sessionThreshold=double(T.sessionThreshold(use));
  common.betaWeibull=betaWeibull;
  common.lapse=lapse;
  common.targetPerformance=targetPerformance;
  common.Xchange=XchangeAll(use,:);
  common.zP=zPAll(use);
  common.zQ=zQAll(use);
  common.nProbeCandidates=nProbeCandidates;

  signedRow=findSignedOffset(signedGainFits,theta);
  change0=signedGainFits.byOffset(signedRow).allSigned.gain(1:4);

  starts6=[ ...
    change0, change0, zeros(4,1); ...
    0.30,    0.60,    0.10; ...
    0.30,    0.60,    0.10];
  lb6=[repmat(gainBounds(1),4,1);0;0];
  ub6=[repmat(gainBounds(2),4,1); ...
    max(0,gainBounds(2));max(0,gainBounds(2))];

  sumFit=fitPoolingModel(sprintf('offset_%g_sum',theta),'sum', ...
    starts6,lb6,ub6,common,opts.PBounds);
  maxStarts=[sumFit.parameter(1:6),starts6];
  hardMaxFit=fitPoolingModel(sprintf('offset_%g_hardMax',theta),'hardMax', ...
    maxStarts,lb6,ub6,common,opts.PBounds);

  pStarts=[ ...
    sumFit.parameter(1:6), hardMaxFit.parameter(1:6), starts6(:,1); ...
    log(1.25),             log(20),                    log(5)];
  lbP=[lb6;log(opts.PBounds(1))];
  ubP=[ub6;log(opts.PBounds(2))];
  pNormFit=fitPoolingModel(sprintf('offset_%g_pNorm',theta),'pNorm', ...
    pStarts,lbP,ubP,common,opts.PBounds);

  localRows=find(use);
  trialNLLSum(localRows)=sumFit.trialNegLogLikelihood;
  trialNLLPNorm(localRows)=pNormFit.trialNegLogLikelihood;
  trialNLLHardMax(localRows)=hardMaxFit.trialNegLogLikelihood;

  byOffset(k).probeDirDeg=theta;
  byOffset(k).nTrials=sumFit.nTrials;
  byOffset(k).nSessions=numel(unique(T.sessionIndex(use)));
  byOffset(k).nProbeCandidates=nProbeCandidates;
  byOffset(k).sum=sumFit;
  byOffset(k).pNorm=pNormFit;
  byOffset(k).hardMax=hardMaxFit;
  byOffset(k).freeSignedNLL=signedGainFits.byOffset(signedRow).allSigned.negLogLikelihood;
  byOffset(k).sessionIndex=unique(T.sessionIndex(use));
end

sharedPNorm=fitSharedPNorm(T,XchangeAll,zPAll,zQAll,byOffset, ...
  betaWeibull,lapse,targetPerformance,gainBounds,opts.PBounds,profilePValues);
trialNLLSharedPNorm=sharedPNorm.trialNegLogLikelihood;

offsetSummary=makeOffsetSummary(byOffset,sharedPNorm);
sessionSummary=makeSessionSummary(T,trialNLLSum,trialNLLPNorm, ...
  trialNLLHardMax,trialNLLSharedPNorm);
overallSummary=makeOverallSummary(byOffset,sharedPNorm,sessionSummary);

F=struct();
F.createdAt=datetime('now');
F.createdBy=mfilename;
F.analysisType='controlled empirical change; opponent-rectified no-change pooling';
F.predictorType='rectStep';
F.preferredNullRatio=opts.PreferredNullRatio;
F.pBounds=opts.PBounds(:)';
F.profilePValues=profilePValues;
F.betaWeibull=betaWeibull;
F.lapse=lapse;
F.targetPerformance=targetPerformance;
F.byOffset=byOffset;
F.sharedPNorm=sharedPNorm;
F.offsetSummary=offsetSummary;
F.sessionSummary=sessionSummary;
F.overallSummary=overallSummary;
F.trialNegLogLikelihood=table( ...
  double(T.sessionIndex),double(T.probeDirDeg),double(T.trialIdx), ...
  trialNLLSum,trialNLLPNorm,trialNLLHardMax,trialNLLSharedPNorm, ...
  'VariableNames',{'sessionIndex','probeDirDeg','trialIdx', ...
  'sumNLL','pNormNLL','hardMaxNLL','sharedPNormNLL'});

if strlength(string(opts.PlotPath))>0
  plotPoolingSummary(F,char(opts.PlotPath));
  F.plotPath=string(opts.PlotPath);
end
end

%% ------------------------------------------------------------------------
function validateInput(D,S,betaWeibull,lapse,targetPerformance,gainBounds)
required={'trialTable','stepFrames','changePrefNoiseByFrameTrial', ...
  'changeProbeEffectiveNoiseByFrameTrial','noChangePrefNoiseByFrameTrial', ...
  'noChangeProbeCandidateNoiseByFrameTrial'};
for k=1:numel(required)
  if ~isfield(D,required{k})
    error('analyzeIDRControlledPooling:MissingField', ...
      'Input sideGainData is missing %s.',required{k});
  end
end
tableRequired={'sessionIndex','trialIdx','probeDirDeg','nYokedProbeStreams', ...
  'stepDeltaCohPC','correct','sessionThreshold'};
missing=setdiff(tableRequired,D.trialTable.Properties.VariableNames);
if ~isempty(missing)
  error('analyzeIDRControlledPooling:MissingTableVariable', ...
    'Trial table is missing: %s.',strjoin(missing,', '));
end
if ~isstruct(S) || ~isfield(S,'byOffset')
  error('analyzeIDRControlledPooling:MissingSignedFits', ...
    'signedGainFits must contain byOffset results.');
end
if ~(isscalar(betaWeibull)&&isfinite(betaWeibull)&&betaWeibull>0)
  error('analyzeIDRControlledPooling:BadBeta','betaWeibull must be positive.');
end
if ~(isscalar(lapse)&&isfinite(lapse)&&lapse>=0&&lapse<0.5)
  error('analyzeIDRControlledPooling:BadLapse','lapse must be in [0,0.5).');
end
if ~(isscalar(targetPerformance)&&targetPerformance>0.5&&targetPerformance<1)
  error('analyzeIDRControlledPooling:BadTarget','targetPerformance must be in (0.5,1).');
end
if ~(isnumeric(gainBounds)&&numel(gainBounds)==2&& ...
    all(isfinite(gainBounds))&&gainBounds(1)<gainBounds(2))
  error('analyzeIDRControlledPooling:BadBounds','gainBounds must be [lower upper].');
end
end

function row=findSignedOffset(S,theta)
values=arrayfun(@(x) double(x.probeDirDeg),S.byOffset);
row=find(values==theta);
if numel(row)~=1
  error('analyzeIDRControlledPooling:SignedOffsetMatch', ...
    'Expected one signed-fit result at offset %g; found %d.',theta,numel(row));
end
end

function x=stepMean(noiseByFrameTrial,stepFrames)
x=mean(double(noiseByFrameTrial(stepFrames,:)),1,'omitnan')';
end
function x=posPart(x)
x=max(double(x),0);
end
function x=negPart(x)
x=min(double(x),0);
end

%% ------------------------------------------------------------------------
function fit=fitPoolingModel(modelName,modelType,starts,lb,ub,common,pBounds)

valid=all(isfinite(common.Xchange),2)&isfinite(common.zP)&isfinite(common.zQ)& ...
  isfinite(common.stepCoh)&isfinite(common.correct)&isfinite(common.sessionThreshold);
C=common;
C.Xchange=common.Xchange(valid,:); C.zP=common.zP(valid); C.zQ=common.zQ(valid);
C.stepCoh=common.stepCoh(valid); C.correct=common.correct(valid);
C.sessionThreshold=common.sessionThreshold(valid);

lb=double(lb(:)); ub=double(ub(:));
if size(starts,1)~=numel(lb)||numel(ub)~=numel(lb)
  error('analyzeIDRControlledPooling:ParameterSize', ...
    'Starts or bounds have the wrong size for %s.',modelName);
end
objective=@(v) poolingNLL(v,modelType,C,pBounds,false);
options=optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',20000,'MaxIterations',5000);

best=nan(numel(lb),1); bestNLL=inf; bestExitflag=nan;
startNLL=nan(size(starts,2),1); exitflags=nan(size(starts,2),1);
for j=1:size(starts,2)
  x0=min(max(starts(:,j),lb),ub);
  [v,nll,exitflag]=fmincon(objective,x0,[],[],[],[],lb,ub,[],options);
  startNLL(j)=nll; exitflags(j)=exitflag;
  if isfinite(nll)&&nll<bestNLL
    best=v(:); bestNLL=nll; bestExitflag=exitflag;
  end
end
if ~isfinite(bestNLL)
  error('analyzeIDRControlledPooling:FitFailure', ...
    'All optimization starts failed for %s.',modelName);
end

H=finiteDifferenceHessian(objective,best);
covariance=nan(numel(best)); SE=nan(numel(best),1);
if all(isfinite(H),'all')
  H=(H+H')/2;
  if all(eig(H)>0)
    covariance=inv(H); covariance=(covariance+covariance')/2;
    variance=diag(covariance); ok=isfinite(variance)&variance>0;
    SE(ok)=sqrt(variance(ok));
  end
end
[~,trialNLL,prediction,poolValue]=poolingNLL(best,modelType,C,pBounds,true);

fit=struct();
fit.model=modelName; fit.modelType=modelType;
fit.fitParameter=best; fit.fitParameterSE=SE;
fit.hessian=H; fit.covariance=covariance;
fit.negLogLikelihood=bestNLL; fit.exitflag=bestExitflag;
fit.startNegLogLikelihoods=startNLL; fit.startExitflags=exitflags;
fit.nParameters=numel(best); fit.nTrials=numel(C.correct);
fit.trialNegLogLikelihood=trialNLL;
fit.predictedCorrectProbability=prediction;
fit.poolValue=poolValue;
fit.parameter=best(1:6);
fit.parameterNames=["gCPpos","gCPneg","gCQpos","gCQneg","aNP","aNQ"];
fit.SE=SE(1:6);
if strcmp(modelType,'pNorm')
  fit.p=exp(best(7));
  fit.logPSE=SE(7);
  fit.pAtLowerBound=abs(fit.p-pBounds(1))<1e-4;
  fit.pAtUpperBound=abs(fit.p-pBounds(2))<1e-3;
else
  fit.p=nan; fit.logPSE=nan; fit.pAtLowerBound=false; fit.pAtUpperBound=false;
end
fit.nEffectiveCohClipped=sum(C.stepCoh + C.Xchange*best(1:4) - poolValue < 0);
if strcmp(modelType,'hardMax')
  uP=best(5)*C.zP; uQ=best(6)*C.zQ;
  fit.probeWinnerFraction=mean(uQ>uP);
  fit.tieFraction=mean(uQ==uP);
else
  fit.probeWinnerFraction=nan; fit.tieFraction=nan;
end
end

%% ------------------------------------------------------------------------
function [nll,trialNLL,pCorrect,poolValue]=poolingNLL(v,modelType,C,pBounds,returnDetails)

changeValue=C.Xchange*v(1:4);
uP=v(5)*C.zP; uQ=v(6)*C.zQ;
switch modelType
  case 'sum'
    poolValue=uP+C.nProbeCandidates*uQ;
  case 'hardMax'
    poolValue=max(uP,uQ);
  case 'pNorm'
    p=min(max(exp(v(7)),pBounds(1)),pBounds(2));
    poolValue=stablePNorm(uP,uQ,C.nProbeCandidates,p);
  otherwise
    error('analyzeIDRControlledPooling:BadModelType','Unknown model type %s.',modelType);
end
effectiveCoh=max(C.stepCoh+changeValue-poolValue,0);
pCorrect=fixedShapeWeibullP(effectiveCoh,C.sessionThreshold, ...
  C.betaWeibull,C.lapse,C.targetPerformance);
pCorrect=min(max(pCorrect,eps),1-eps);
trialNLL=-(C.correct.*log(pCorrect)+(1-C.correct).*log(1-pCorrect));
nll=sum(trialNLL);
if ~returnDetails
  trialNLL=[]; pCorrect=[];
end
end

function y=stablePNorm(uP,uQ,nProbeCandidates,p)
m=max(uP,uQ); y=zeros(size(m)); use=m>0;
a=zeros(size(m)); b=zeros(size(m));
a(use)=uP(use)./m(use); b(use)=uQ(use)./m(use);
y(use)=m(use).*(a(use).^p+nProbeCandidates*b(use).^p).^(1/p);
end

function p=fixedShapeWeibullP(coh,threshold,betaWeibull,lapse,targetPerformance)
ratio=(1-lapse-targetPerformance)/(0.5-lapse);
alpha=threshold./((-log(ratio))^(1/betaWeibull));
p=1-lapse-(0.5-lapse).*exp(-(max(double(coh),0)./alpha).^betaWeibull);
end

function H=finiteDifferenceHessian(objective,x)
x=x(:); n=numel(x); H=nan(n); h=3e-4.*max(1,abs(x)); f0=objective(x);
for i=1:n
  ei=zeros(n,1); ei(i)=h(i); fp=objective(x+ei); fm=objective(x-ei);
  H(i,i)=(fp-2*f0+fm)/h(i)^2;
  for j=i+1:n
    ej=zeros(n,1); ej(j)=h(j);
    fpp=objective(x+ei+ej); fpm=objective(x+ei-ej);
    fmp=objective(x-ei+ej); fmm=objective(x-ei-ej);
    H(i,j)=(fpp-fpm-fmp+fmm)/(4*h(i)*h(j)); H(j,i)=H(i,j);
  end
end
end

%% ------------------------------------------------------------------------
function fit=fitSharedPNorm(T,XchangeAll,zPAll,zQAll,R, ...
  betaWeibull,lapse,targetPerformance,gainBounds,pBounds,profilePValues)
% Fit one p across offsets while retaining six offset-specific parameters.
% A full finite-difference Hessian is intentionally omitted: it is costly
% for this 43-parameter descriptive comparison and is not needed for AIC.

nOffsets=numel(R); nPerOffset=6; nParameters=nPerOffset*nOffsets+1;
C=cell(nOffsets,1); validRows=cell(nOffsets,1);
for k=1:nOffsets
  use=double(T.probeDirDeg)==R(k).probeDirDeg;
  valid=use & all(isfinite(XchangeAll),2) & isfinite(zPAll) & ...
    isfinite(zQAll) & isfinite(double(T.stepDeltaCohPC)) & ...
    isfinite(double(T.correct)) & isfinite(double(T.sessionThreshold));
  validRows{k}=find(valid);
  C{k}=struct( ...
    'stepCoh',double(T.stepDeltaCohPC(valid)), ...
    'correct',double(T.correct(valid)), ...
    'sessionThreshold',double(T.sessionThreshold(valid)), ...
    'betaWeibull',betaWeibull,'lapse',lapse, ...
    'targetPerformance',targetPerformance, ...
    'Xchange',XchangeAll(valid,:), ...
    'zP',zPAll(valid),'zQ',zQAll(valid), ...
    'nProbeCandidates',R(k).nProbeCandidates);
end

lb6=[repmat(gainBounds(1),4,1);0;0];
ub6=[repmat(gainBounds(2),4,1); ...
  max(0,gainBounds(2));max(0,gainBounds(2))];
lb=[repmat(lb6,nOffsets,1);log(pBounds(1))];
ub=[repmat(ub6,nOffsets,1);log(pBounds(2))];

starts=nan(nParameters,4);
for k=1:nOffsets
  ii=(k-1)*nPerOffset+(1:nPerOffset);
  starts(ii,1)=R(k).pNorm.parameter;
  starts(ii,2)=R(k).sum.parameter;
  starts(ii,3)=R(k).hardMax.parameter;
  starts(ii,4)=R(k).hardMax.parameter;
end
separateP=arrayfun(@(x)x.pNorm.p,R);
starts(end,:)=[median(log(separateP(isfinite(separateP)))), ...
  log(1.01),log(20),log(80)];
starts=min(max(starts,lb),ub);

objective=@(v) sharedPNormNLL(v,C,pBounds,false);
options=optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',50000,'MaxIterations',5000);
best=nan(nParameters,1); bestNLL=inf; bestExitflag=nan;
startNLL=nan(size(starts,2),1); exitflags=nan(size(starts,2),1);
for j=1:size(starts,2)
  [v,nll,exitflag]=fmincon(objective,starts(:,j),[],[],[],[], ...
    lb,ub,[],options);
  startNLL(j)=nll; exitflags(j)=exitflag;
  if isfinite(nll) && nll<bestNLL
    best=v(:); bestNLL=nll; bestExitflag=exitflag;
  end
end
if ~isfinite(bestNLL)
  error('analyzeIDRControlledPooling:SharedFitFailure', ...
    'All optimization starts failed for the shared-p model.');
end

[~,trialByOffset,offsetNLL]=sharedPNormNLL(best,C,pBounds,true);
trialNLL=nan(height(T),1);
for k=1:nOffsets
  trialNLL(validRows{k})=trialByOffset{k};
end
parameterByOffset=reshape(best(1:nPerOffset*nOffsets),nPerOffset,nOffsets)';
profile=fitSharedPProfile(C,R,parameterByOffset,bestNLL, ...
  lb6,ub6,pBounds,profilePValues);
fit=struct();
fit.model='shared_pNorm'; fit.modelType='pNorm';
fit.fitParameter=best; fit.negLogLikelihood=bestNLL;
fit.exitflag=bestExitflag; fit.startNegLogLikelihoods=startNLL;
fit.startExitflags=exitflags; fit.nParameters=nParameters;
fit.nTrials=sum(cellfun(@numel,validRows));
fit.p=exp(best(end)); fit.logP=best(end);
fit.pAtLowerBound=abs(fit.p-pBounds(1))<1e-4;
fit.pAtUpperBound=abs(fit.p-pBounds(2))<1e-3;
fit.parameterByOffset=parameterByOffset;
fit.parameterNames=["gCPpos","gCPneg","gCQpos","gCQneg","aNP","aNQ"];
fit.probeDirDeg=reshape(arrayfun(@(x)x.probeDirDeg,R),[],1);
fit.parameterTable=array2table(parameterByOffset, ...
  'VariableNames',cellstr(fit.parameterNames));
fit.parameterTable=addvars(fit.parameterTable,fit.probeDirDeg, ...
  'Before',1,'NewVariableNames','probeDirDeg');
fit.offsetNegLogLikelihood=offsetNLL;
fit.trialNegLogLikelihood=trialNLL;
fit.profile=profile;
end

function [nll,trialByOffset,offsetNLL]=sharedPNormNLL(v,C,pBounds,returnDetails)
nOffsets=numel(C); nPerOffset=6; sharedLogP=v(end);
offsetNLL=nan(nOffsets,1); trialByOffset=cell(nOffsets,1);
for k=1:nOffsets
  ii=(k-1)*nPerOffset+(1:nPerOffset);
  local=[v(ii);sharedLogP];
  [offsetNLL(k),trialByOffset{k}]=poolingNLL( ...
    local,'pNorm',C{k},pBounds,returnDetails);
end
nll=sum(offsetNLL);
if ~returnDetails
  trialByOffset=[]; offsetNLL=[];
end
end

function profile=fitSharedPProfile(C,R,sharedParameterByOffset, ...
  unconstrainedNLL,lb6,ub6,pBounds,pValues)
% At fixed p, offsets are separable. Reoptimize all six nuisance parameters
% independently by offset, using the preceding p and shared fit as starts.
nP=numel(pValues); nOffsets=numel(C); nParametersPerOffset=6;
parameterByPAndOffset=nan(nP,nOffsets,nParametersPerOffset);
offsetNLL=nan(nP,nOffsets); exitflag=nan(nP,nOffsets);
startNLL=nan(nP,nOffsets,2); startExitflag=nan(nP,nOffsets,2);
warm=nan(nParametersPerOffset,nOffsets);
for k=1:nOffsets
  warm(:,k)=R(k).sum.parameter;
end
options=optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',20000,'MaxIterations',5000);

for i=1:nP
  fprintf('  Shared-p profile: p = %g (%d/%d)\n',pValues(i),i,nP);
  fixedLogP=log(pValues(i));
  for k=1:nOffsets
    starts=[warm(:,k),sharedParameterByOffset(k,:)'];
    objective=@(w) poolingNLL([w(:);fixedLogP], ...
      'pNorm',C{k},pBounds,false);
    best=nan(nParametersPerOffset,1); bestNLL=inf; bestExit=nan;
    for j=1:size(starts,2)
      x0=min(max(starts(:,j),lb6),ub6);
      [v,nll,thisExit]=fmincon(objective,x0,[],[],[],[], ...
        lb6,ub6,[],options);
      startNLL(i,k,j)=nll; startExitflag(i,k,j)=thisExit;
      if isfinite(nll) && nll<bestNLL
        best=v(:); bestNLL=nll; bestExit=thisExit;
      end
    end
    if ~isfinite(bestNLL)
      error('analyzeIDRControlledPooling:ProfileFitFailure', ...
        'All profile starts failed at p=%g, offset=%g.', ...
        pValues(i),R(k).probeDirDeg);
    end
    warm(:,k)=best;
    parameterByPAndOffset(i,k,:)=reshape(best,1,1,[]);
    offsetNLL(i,k)=bestNLL; exitflag(i,k)=bestExit;
  end
end

negLogLikelihood=sum(offsetNLL,2);
deltaNLL=max(negLogLikelihood-unconstrainedNLL,0);
relativeLikelihood=exp(-deltaNLL);
cutoffOneSided95=1.352772;  % 0.5*chi2inv(0.90,1)
cutoffTwoSided95=1.920729;  % 0.5*chi2inv(0.95,1)
withinOneSided95=deltaNLL<=cutoffOneSided95;
withinTwoSided95=deltaNLL<=cutoffTwoSided95;
[lowerOneSided95,lowerOneTruncated]=profileLowerBound( ...
  pValues,deltaNLL,cutoffOneSided95);
[lowerTwoSided95,lowerTwoTruncated]=profileLowerBound( ...
  pValues,deltaNLL,cutoffTwoSided95);

profileTable=table(pValues,negLogLikelihood,deltaNLL,relativeLikelihood, ...
  withinOneSided95,withinTwoSided95, ...
  'VariableNames',{'p','negLogLikelihood','deltaNLL','relativeLikelihood', ...
  'withinOneSided95','withinTwoSided95'});
for k=1:nOffsets
  variableName=sprintf('offset%gNLL',R(k).probeDirDeg);
  profileTable.(variableName)=offsetNLL(:,k);
end

profile=struct();
profile.table=profileTable;
profile.pValues=pValues;
profile.negLogLikelihood=negLogLikelihood;
profile.deltaNLL=deltaNLL;
profile.relativeLikelihood=relativeLikelihood;
profile.offsetNegLogLikelihood=offsetNLL;
profile.parameterByPAndOffset=parameterByPAndOffset;
profile.exitflag=exitflag;
profile.startNegLogLikelihoods=startNLL;
profile.startExitflags=startExitflag;
profile.nFittedParameters=numel(R)*nParametersPerOffset;
profile.cutoffOneSided95=cutoffOneSided95;
profile.cutoffTwoSided95=cutoffTwoSided95;
profile.lowerOneSided95=lowerOneSided95;
profile.lowerTwoSided95=lowerTwoSided95;
profile.lowerOneSided95GridTruncated=lowerOneTruncated;
profile.lowerTwoSided95GridTruncated=lowerTwoTruncated;
sumEndpoint=find(abs(pValues-1)<1e-12,1);
maxEndpoint=find(abs(pValues-pBounds(2))<1e-12,1);
profile.sumEndpointDifference=nan;
profile.maxEndpointDifference=nan;
if ~isempty(sumEndpoint)
  profile.sumEndpointDifference=negLogLikelihood(sumEndpoint)- ...
    sum(arrayfun(@(x)x.sum.negLogLikelihood,R));
end
if ~isempty(maxEndpoint)
  profile.maxEndpointDifference=negLogLikelihood(maxEndpoint)- ...
    sum(arrayfun(@(x)x.hardMax.negLogLikelihood,R));
end
end

function [lowerBound,gridTruncated]=profileLowerBound(pValues,deltaNLL,cutoff)
firstInside=find(deltaNLL<=cutoff,1,'first');
gridTruncated=false;
if isempty(firstInside)
  lowerBound=nan;
elseif firstInside==1
  lowerBound=pValues(1); gridTruncated=true;
else
  i0=firstInside-1; i1=firstInside;
  d0=deltaNLL(i0); d1=deltaNLL(i1);
  fraction=(cutoff-d0)/(d1-d0);
  lowerBound=exp(log(pValues(i0))+fraction* ...
    (log(pValues(i1))-log(pValues(i0))));
end
end

%% ------------------------------------------------------------------------
function S=makeOffsetSummary(R,sharedPNorm)
n=numel(R);
probeDirDeg=nan(n,1); nTrials=zeros(n,1); nSessions=zeros(n,1);
nProbeCandidates=zeros(n,1); freeSignedNLL=nan(n,1);
sumNLL=nan(n,1); pNormNLL=nan(n,1); hardMaxNLL=nan(n,1);
sharedPNormNLL=sharedPNorm.offsetNegLogLikelihood;
deltaNLLSumMinusMax=nan(n,1); deltaNLLSumMinusPNorm=nan(n,1);
deltaNLLMaxMinusPNorm=nan(n,1); fittedP=nan(n,1);
deltaNLLSumMinusShared=nan(n,1); deltaNLLMaxMinusShared=nan(n,1);
sumAIC=nan(n,1); pNormAIC=nan(n,1); hardMaxAIC=nan(n,1);
maxProbeWinnerFraction=nan(n,1); maxTieFraction=nan(n,1);
sumANP=nan(n,1); sumANQ=nan(n,1); maxANP=nan(n,1); maxANQ=nan(n,1);
sharedP=repmat(sharedPNorm.p,n,1);
sharedGCPpos=nan(n,1); sharedGCPneg=nan(n,1);
sharedGCQpos=nan(n,1); sharedGCQneg=nan(n,1);
sharedANP=nan(n,1); sharedANQ=nan(n,1);
for k=1:n
  probeDirDeg(k)=R(k).probeDirDeg; nTrials(k)=R(k).nTrials;
  nSessions(k)=R(k).nSessions; nProbeCandidates(k)=R(k).nProbeCandidates;
  freeSignedNLL(k)=R(k).freeSignedNLL;
  sumNLL(k)=R(k).sum.negLogLikelihood;
  pNormNLL(k)=R(k).pNorm.negLogLikelihood;
  hardMaxNLL(k)=R(k).hardMax.negLogLikelihood;
  deltaNLLSumMinusMax(k)=sumNLL(k)-hardMaxNLL(k);
  deltaNLLSumMinusPNorm(k)=sumNLL(k)-pNormNLL(k);
  deltaNLLMaxMinusPNorm(k)=hardMaxNLL(k)-pNormNLL(k);
  deltaNLLSumMinusShared(k)=sumNLL(k)-sharedPNormNLL(k);
  deltaNLLMaxMinusShared(k)=hardMaxNLL(k)-sharedPNormNLL(k);
  fittedP(k)=R(k).pNorm.p;
  sumAIC(k)=2*sumNLL(k)+2*R(k).sum.nParameters;
  pNormAIC(k)=2*pNormNLL(k)+2*R(k).pNorm.nParameters;
  hardMaxAIC(k)=2*hardMaxNLL(k)+2*R(k).hardMax.nParameters;
  maxProbeWinnerFraction(k)=R(k).hardMax.probeWinnerFraction;
  maxTieFraction(k)=R(k).hardMax.tieFraction;
  sumANP(k)=R(k).sum.parameter(5); sumANQ(k)=R(k).sum.parameter(6);
  maxANP(k)=R(k).hardMax.parameter(5); maxANQ(k)=R(k).hardMax.parameter(6);
  sharedGCPpos(k)=sharedPNorm.parameterByOffset(k,1);
  sharedGCPneg(k)=sharedPNorm.parameterByOffset(k,2);
  sharedGCQpos(k)=sharedPNorm.parameterByOffset(k,3);
  sharedGCQneg(k)=sharedPNorm.parameterByOffset(k,4);
  sharedANP(k)=sharedPNorm.parameterByOffset(k,5);
  sharedANQ(k)=sharedPNorm.parameterByOffset(k,6);
end
S=table(probeDirDeg,nTrials,nSessions,nProbeCandidates,freeSignedNLL, ...
  sumNLL,pNormNLL,hardMaxNLL,sharedPNormNLL,deltaNLLSumMinusMax, ...
  deltaNLLSumMinusPNorm,deltaNLLMaxMinusPNorm,fittedP, ...
  deltaNLLSumMinusShared,deltaNLLMaxMinusShared, ...
  sumAIC,pNormAIC,hardMaxAIC,maxProbeWinnerFraction,maxTieFraction, ...
  sumANP,sumANQ,maxANP,maxANQ,sharedP,sharedGCPpos,sharedGCPneg, ...
  sharedGCQpos,sharedGCQneg,sharedANP,sharedANQ);
end

function S=makeSessionSummary(T,nllSum,nllPNorm,nllMax,nllShared)
[G,sessionIndex,sessionID]=findgroups(double(T.sessionIndex),string(T.sessionID));
nTrials=splitapply(@numel,nllSum,G);
sumNLL=splitapply(@sum,nllSum,G);
pNormNLL=splitapply(@sum,nllPNorm,G);
hardMaxNLL=splitapply(@sum,nllMax,G);
sharedPNormNLL=splitapply(@sum,nllShared,G);
deltaNLLSumMinusMax=sumNLL-hardMaxNLL;
deltaNLLSumMinusPNorm=sumNLL-pNormNLL;
deltaNLLSumMinusShared=sumNLL-sharedPNormNLL;
deltaNLLMaxMinusShared=hardMaxNLL-sharedPNormNLL;
S=table(sessionIndex,sessionID,nTrials,sumNLL,pNormNLL,hardMaxNLL, ...
  sharedPNormNLL,deltaNLLSumMinusMax,deltaNLLSumMinusPNorm, ...
  deltaNLLSumMinusShared,deltaNLLMaxMinusShared);
end

function S=makeOverallSummary(R,sharedPNorm,sessionSummary)
sumNLL=sum(arrayfun(@(x)x.sum.negLogLikelihood,R));
pNormNLL=sum(arrayfun(@(x)x.pNorm.negLogLikelihood,R));
hardMaxNLL=sum(arrayfun(@(x)x.hardMax.negLogLikelihood,R));
sharedPNormNLL=sharedPNorm.negLogLikelihood;
nParametersSum=sum(arrayfun(@(x)x.sum.nParameters,R));
nParametersPNorm=sum(arrayfun(@(x)x.pNorm.nParameters,R));
nParametersHardMax=sum(arrayfun(@(x)x.hardMax.nParameters,R));
nParametersSharedPNorm=sharedPNorm.nParameters;
deltaNLLSumMinusMax=sumNLL-hardMaxNLL;
deltaNLLSumMinusPNorm=sumNLL-pNormNLL;
deltaNLLMaxMinusPNorm=hardMaxNLL-pNormNLL;
deltaNLLSumMinusShared=sumNLL-sharedPNormNLL;
deltaNLLMaxMinusShared=hardMaxNLL-sharedPNormNLL;
sumAIC=2*sumNLL+2*nParametersSum;
pNormAIC=2*pNormNLL+2*nParametersPNorm;
hardMaxAIC=2*hardMaxNLL+2*nParametersHardMax;
sharedPNormAIC=2*sharedPNormNLL+2*nParametersSharedPNorm;
nSessions=height(sessionSummary);
nSessionsFavorMax=sum(sessionSummary.deltaNLLSumMinusMax>0);
nSessionsFavorSum=sum(sessionSummary.deltaNLLSumMinusMax<0);
nSessionsFavorSharedVsMax=sum(sessionSummary.deltaNLLMaxMinusShared>0);
nSessionsFavorMaxVsShared=sum(sessionSummary.deltaNLLMaxMinusShared<0);
sharedP=sharedPNorm.p;
profileLowerOneSided95=sharedPNorm.profile.lowerOneSided95;
profileLowerTwoSided95=sharedPNorm.profile.lowerTwoSided95;
S=table(sumNLL,pNormNLL,hardMaxNLL,sharedPNormNLL, ...
  nParametersSum,nParametersPNorm,nParametersHardMax,nParametersSharedPNorm, ...
  deltaNLLSumMinusMax,deltaNLLSumMinusPNorm,deltaNLLMaxMinusPNorm, ...
  deltaNLLSumMinusShared,deltaNLLMaxMinusShared, ...
  sumAIC,pNormAIC,hardMaxAIC,sharedPNormAIC,sharedP, ...
  profileLowerOneSided95,profileLowerTwoSided95,nSessions, ...
  nSessionsFavorMax,nSessionsFavorSum,nSessionsFavorSharedVsMax, ...
  nSessionsFavorMaxVsShared);
end

%% ------------------------------------------------------------------------
function plotPoolingSummary(F,pdfPath)
S=F.offsetSummary; P=F.sharedPNorm.profile;
fig=figure('Color','w','Visible','off','Units','inches','Position',[1 1 10 10]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('IDR controlled no-change pooling; preferred:null = %.2f', ...
  F.preferredNullRatio),'FontWeight','bold');
ax=nexttile(tl); hold(ax,'on');
plot(ax,S.probeDirDeg,S.deltaNLLSumMinusMax,'o-','LineWidth',1.3);
yline(ax,0,'k:'); ylabel(ax,'NLL(sum) - NLL(max)');
title(ax,'Positive values favor hard max'); grid(ax,'on'); box(ax,'off');
ax=nexttile(tl); plot(ax,S.probeDirDeg,S.fittedP,'o-','LineWidth',1.3);
yline(ax,F.sharedPNorm.p,'--','Shared p','LineWidth',1.1);
yline(ax,1,'k:'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Fitted p');
title(ax,'Opponent-rectified p-norm'); grid(ax,'on'); box(ax,'off');
ax=nexttile(tl); hold(ax,'on');
semilogx(ax,P.pValues,P.deltaNLL,'o-','LineWidth',1.3);
yline(ax,P.cutoffOneSided95,'--','One-sided 95% lower bound');
yline(ax,P.cutoffTwoSided95,':','Two-sided 95% interval');
xlabel(ax,'Fixed shared p'); ylabel(ax,'Profile \DeltaNLL');
title(ax,'Shared-p profile likelihood'); grid(ax,'on'); box(ax,'off');
exportgraphics(fig,pdfPath,'ContentType','vector'); close(fig);
end
