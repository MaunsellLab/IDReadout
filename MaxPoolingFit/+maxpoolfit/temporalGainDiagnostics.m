function diagnostics = temporalGainDiagnostics(cache,prediction,varargin)
%TEMPORALGAINDIAGNOSTICS Five-bin changed/no-change RDK noise regression.
%
% The cached 120-Hz step traces are averaged into simultaneous 50-ms
% predictors. Changed- and unchanged-RDK drift and probe coefficients are
% fitted jointly with session fixed effects. IDR uses one copy of its yoked
% probe trace; IDQ/IDS use the sum of their two independent probe traces.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('Domain',"",@(x)isstring(x)||ischar(x));
p.addParameter('Animal',"",@(x)isstring(x)||ischar(x));
p.addParameter('ProbeOffsetDeg',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
p.addParameter('StepMagnitudePC',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
p.addParameter('BinFrames',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.addParameter('FrameRateHz',120,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinimumTrials',500,@(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('MinimumSessions',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('NBootstrap',200,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&mod(x,1)==0);
p.addParameter('RandomSeed',271828,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.parse(varargin{:});

maxpoolfit.validateCache(cache);
assert(isfield(cache,'raw')&&isfield(cache.raw,'stepStreamCohPC'), ...
  'The cache lacks raw framewise step-stream coherence.');
T=cache.trials; mask=trialMask(T,p.Results.TrialMask);
domain=string(p.Results.Domain); animal=string(p.Results.Animal);
if strlength(domain)>0, mask=mask & string(T.domain)==domain; end
if strlength(animal)>0, mask=mask & string(T.animal)==animal; end
if ~isempty(p.Results.ProbeOffsetDeg)
  mask=mask & ismember(T.probeOffsetDeg,double(p.Results.ProbeOffsetDeg));
end
if ~isempty(p.Results.StepMagnitudePC)
  mask=mask & ismember(T.stepMagnitudePC,double(p.Results.StepMagnitudePC));
end
idx=find(mask);
assert(numel(idx)>=p.Results.MinimumTrials, ...
  'Only %d selected trials; at least %d are required.',numel(idx),p.Results.MinimumTrials);
assert(numel(prediction.pRDK1)==height(T), ...
  'prediction.pRDK1 must contain one value per cached trial.');
assert(all(isfinite(prediction.pRDK1(idx))), ...
  'Prediction does not cover every selected trial.');

D=framePredictors(cache,idx,p.Results.BinFrames,p.Results.FrameRateHz);
[~,sessions]=findgroups(D.sessionID);
assert(numel(sessions)>=p.Results.MinimumSessions, ...
  'Only %d selected sessions; at least %d are required.',numel(sessions),p.Results.MinimumSessions);

observedCorrect=double(T.choiceRDK(idx)==T.changeRDK(idx));
change1=double(T.changeRDK(idx)==1);
predictedCorrect=change1.*double(prediction.pRDK1(idx))+ ...
  (1-change1).*(1-double(prediction.pRDK1(idx)));
outcomes=[observedCorrect predictedCorrect];
sources=["observed" "primary domain noise"];

old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(p.Results.RandomSeed,'twister');
bootstrapSamples=makeBootstrapSamples(D.sessionID,p.Results.NBootstrap);
profile=table(); shapeTests=table(); fits=cell(1,numel(sources));
for s=1:numel(sources)
  [b,V,fitInfo,columns,predictorAudit]=fitTemporal(D,outcomes(:,s));
  fits{s}=fitInfo;
  [bootChanged,bootUnchanged]=bootstrapProfiles(D,outcomes(:,s),bootstrapSamples,D.nBins);
  [sourceProfile,sourceTests]=summarizeSource( ...
    sources(s),b,V,columns,bootChanged,bootUnchanged,D, ...
    mean(outcomes(:,s)),numel(sessions),predictorAudit);
  profile=[profile;sourceProfile]; %#ok<AGROW>
  shapeTests=[shapeTests;sourceTests]; %#ok<AGROW>
end

selection=struct('domain',unique(string(T.domain(idx)),'stable'), ...
  'animal',unique(string(T.animal(idx)),'stable'), ...
  'probeOffsetDeg',unique(T.probeOffsetDeg(idx),'sorted'), ...
  'stepMagnitudePC',unique(T.stepMagnitudePC(idx),'sorted'), ...
  'nTrials',numel(idx),'nSessions',numel(sessions));
diagnostics=struct('profile',profile,'shapeTests',shapeTests, ...
  'selection',selection,'trialIdx',idx,'nTrials',numel(idx), ...
  'nSessions',numel(sessions),'nFrames',D.nFrames,'nBins',D.nBins, ...
  'binFrames',p.Results.BinFrames,'frameRateHz',p.Results.FrameRateHz, ...
  'binDurationMS',1000*p.Results.BinFrames/p.Results.FrameRateHz, ...
  'nBootstrap',p.Results.NBootstrap,'randomSeed',p.Results.RandomSeed, ...
  'sources',sources,'fitInfo',{fits}, ...
  'regressionDefinition',[ ...
    "joint five-bin fractional-binomial regression; session fixed effects; " ...
    "change-side and within-session step controls; within-session centered noise predictors"], ...
  'probePredictorContract',[ ...
    "IDR: one yoked probe trace; IDQ/IDS: sum of two independent probe traces"], ...
  'createdAt',datetime('now'));
end

function mask=trialMask(T,input)
if isempty(input)
  mask=T.valid & T.hasStepNoise & T.isIncrement;
elseif islogical(input)
  assert(numel(input)==height(T),'Logical TrialMask must contain one value per cached trial.');
  mask=input(:);
else
  mask=false(height(T),1); mask(input)=true;
end
end

function D=framePredictors(cache,idx,binFrames,frameRateHz)
T=cache.trials(idx,:); n=numel(idx); frameCounts=zeros(n,2);
for k=1:n
  for j=1:2
    C=cache.raw.stepStreamCohPC{idx(k),j};
    assert(isnumeric(C)&&ismatrix(C)&&all(isfinite(C),'all'), ...
      'Framewise coherence must be a finite numeric matrix.');
    frameCounts(k,j)=size(C,1);
  end
end
assert(numel(unique(frameCounts(:)))==1,'Selected trials do not share one step-frame count.');
nFrames=frameCounts(1);
if mod(nFrames,binFrames)~=0
  error('maxpoolfit:temporalGainDiagnostics:NonDivisibleFrameCount', ...
    '%d step frames cannot be divided into bins of %d frames.',nFrames,binFrames);
end
nBins=nFrames/binFrames; assert(nBins>=2,'Temporal diagnostics require at least two bins.');
X=nan(n,4*nBins);
for k=1:n
  c=double(T.changeRDK(k)); u=3-c;
  changed=double(cache.raw.stepStreamCohPC{idx(k),c});
  unchanged=double(cache.raw.stepStreamCohPC{idx(k),u});
  assert(size(changed,2)==numel(cache.streamDirDeg{idx(k),c}) && ...
    size(unchanged,2)==numel(cache.streamDirDeg{idx(k),u}), ...
    'Framewise stream count does not match stream directions.');
  assert(size(changed,2)>=2&&size(unchanged,2)>=2, ...
    'Every selected RDK must contain drift and probe streams.');
  for b=1:nBins
    rows=(b-1)*binFrames+(1:binFrames);
    cm=mean(changed(rows,:),1); um=mean(unchanged(rows,:),1);
    X(k,b)=cm(1)-double(T.stepEndpointCohPC(k));
    X(k,nBins+b)=um(1)-double(T.baseCohPC(k));
    if string(T.domain(k))=="IDR"
      X(k,2*nBins+b)=cm(2); X(k,3*nBins+b)=um(2);
    else
      X(k,2*nBins+b)=sum(cm(2:end)); X(k,3*nBins+b)=sum(um(2:end));
    end
  end
end
D=struct('X',X,'sessionID',string(T.sessionID), ...
  'changeRDK',double(T.changeRDK),'stepMagnitudePC',double(T.stepMagnitudePC), ...
  'nFrames',nFrames,'nBins',nBins,'binFrames',binFrames,'frameRateHz',frameRateHz);
end

function samples=makeBootstrapSamples(sessionID,nBootstrap)
[group,sessions]=findgroups(string(sessionID)); nSessions=numel(sessions);
samples=cell(nBootstrap,2);
for r=1:nBootstrap
  draws=randi(nSessions,nSessions,1); rows=[]; labels=strings(0,1);
  for k=1:nSessions
    q=find(group==draws(k)); rows=[rows;q]; %#ok<AGROW>
    labels=[labels;repmat("bootstrap"+r+"_session"+k,numel(q),1)]; %#ok<AGROW>
  end
  samples{r,1}=rows; samples{r,2}=labels;
end
end

function [b,V,info,columns,audit]=fitTemporal(D,y,varargin)
if nargin>=3, sessionID=varargin{1}; else, sessionID=D.sessionID; end
[X,columns,cluster,audit]=designMatrix(D,sessionID);
[b,V,info]=fitBinomial(X,y,cluster);
end

function [X,columns,sessionGroup,audit]=designMatrix(D,sessionID)
n=size(D.X,1); [sessionGroup,sessions]=findgroups(string(sessionID));
nSessions=numel(sessions);
sessionIntercept=sparse((1:n)',sessionGroup,1,n,nSessions);
sideSign=2.*double(D.changeRDK==1)-1;
step=centerWithinGroup(double(D.stepMagnitudePC),sessionGroup);
noise=double(D.X);
for k=1:size(noise,2), noise(:,k)=centerWithinGroup(noise(:,k),sessionGroup); end
nuisance=sparse(sideSign);
if any(abs(step)>1e-10), nuisance=[nuisance sparse(step)]; end
X=[sessionIntercept nuisance sparse(noise)];
columns=(size(X,2)-size(noise,2)+1):size(X,2);
scale=std(noise,0,1); q=scale>1e-12;
if sum(q)>=2
  R=corrcoef(noise(:,q)); predictorConditionNumber=cond(R);
  maximumAbsoluteCorrelation=max(abs(R(~eye(size(R)))));
else
  predictorConditionNumber=nan; maximumAbsoluteCorrelation=nan;
end
audit=struct('nSessions',nSessions,'predictorConditionNumber',predictorConditionNumber, ...
  'maximumAbsolutePredictorCorrelation',maximumAbsoluteCorrelation);
end

function centered=centerWithinGroup(x,group)
means=splitapply(@mean,x,group); centered=x-means(group);
end

function [changed,unchanged]=bootstrapProfiles(D,y,samples,nBins)
nBootstrap=size(samples,1);
changed=nan(nBootstrap,2,nBins); unchanged=nan(nBootstrap,2,nBins);
for r=1:nBootstrap
  rows=samples{r,1}; B=subsetData(D,rows);
  [b,~,~,columns]=fitTemporal(B,y(rows),samples{r,2});
  coefficient=b(columns);
  for component=1:2
    cidx=(component-1)*2*nBins+(1:nBins); uidx=cidx+nBins;
    c=coefficient(cidx); u=-coefficient(uidx);
    if abs(sum(c))>1e-10
      changed(r,component,:)=reshape(c./sum(c),1,1,nBins);
    end
    if abs(sum(u))>1e-10
      unchanged(r,component,:)=reshape(u./sum(u),1,1,nBins);
    end
  end
end
end

function B=subsetData(D,rows)
B=D; B.X=D.X(rows,:); B.sessionID=D.sessionID(rows);
B.changeRDK=D.changeRDK(rows); B.stepMagnitudePC=D.stepMagnitudePC(rows);
end

function [P,S]=summarizeSource(source,b,V,columns,bootC,bootU,D,accuracy,nSessions,audit)
nBins=D.nBins; coefficient=b(columns); covariance=V(columns,columns);
P=table(); S=table(); components=["drift" "probe"];
binDurationMS=1000*D.binFrames/D.frameRateHz;
for component=1:2
  cidx=(component-1)*2*nBins+(1:nBins); uidx=cidx+nBins;
  c=coefficient(cidx); u=coefficient(uidx); cs=c; us=-u;
  cn=normalizeProfile(cs); un=normalizeProfile(us);
  [wald,df,pValue]=shapeWald(c,u,covariance([cidx uidx],[cidx uidx]));
  changedTotal=sum(cs); unchangedTotal=sum(us);
  if abs(changedTotal)>1e-10, ratio=unchangedTotal/changedTotal; else, ratio=nan; end
  sourceRow=table(source,components(component),size(D.X,1),nSessions,accuracy, ...
    changedTotal,unchangedTotal,ratio,wald,df,pValue, ...
    audit.predictorConditionNumber,audit.maximumAbsolutePredictorCorrelation, ...
    'VariableNames',{'source','component','nTrials','nSessions','accuracy', ...
    'changedWholeStepCoefficientPerPC','unchangedWholeStepCoefficientPerPC', ...
    'unchangedOverChanged','shapeWaldStatistic','shapeDF','shapePValue', ...
    'predictorConditionNumber','maximumAbsolutePredictorCorrelation'});
  S=[S;sourceRow]; %#ok<AGROW>
  for status=1:2
    if status==1
      rdkStatus="changed"; raw=cs; rawSE=sqrt(max(diag(covariance(cidx,cidx)),0));
      normalized=cn; boot=reshape(bootC(:,component,:),size(bootC,1),nBins);
    else
      rdkStatus="unchanged"; raw=us; rawSE=sqrt(max(diag(covariance(uidx,uidx)),0));
      normalized=un; boot=reshape(bootU(:,component,:),size(bootU,1),nBins);
    end
    for q=1:nBins
      timeBin=q; timeStartMS=(q-1)*binDurationMS; timeEndMS=q*binDurationMS;
      [lower,upper,nValid]=percentileInterval(boot(:,q));
      row=table(source,components(component),rdkStatus,timeBin,timeStartMS,timeEndMS, ...
        raw(q),rawSE(q),normalized(q),lower,upper,nValid, ...
        'VariableNames',{'source','component','rdkStatus','timeBin', ...
        'timeStartMS','timeEndMS','signedCoefficientPerPC','coefficientSE', ...
        'normalizedCoefficient','bootstrapLower95','bootstrapUpper95','nValidBootstrap'});
      P=[P;row]; %#ok<AGROW>
    end
  end
end
end

function value=normalizeProfile(value)
total=sum(value);
if abs(total)>1e-10, value=value./total; else, value(:)=nan; end
end

function [statistic,df,pValue]=shapeWald(c,u,V)
n=numel(c); sc=sum(c); su=sum(u); statistic=nan; df=nan; pValue=nan;
if abs(sc)<=1e-10||abs(su)<=1e-10||any(~isfinite(V),'all'), return; end
d=c./sc-u./su;
Jc=(eye(n).*sc-c*ones(1,n))./(sc.^2);
Ju=-(eye(n).*su-u*ones(1,n))./(su.^2);
J=[Jc Ju]; W=J*V*J';
df=rank(W,1e-9*max(1,norm(W,2)));
if df<1, statistic=nan; df=nan; return; end
statistic=max(0,d'*pinv(W)*d); pValue=gammainc(statistic/2,df/2,'upper');
end

function [lower,upper,n]=percentileInterval(x)
x=sort(x(isfinite(x))); n=numel(x);
if n<2, lower=nan; upper=nan; return; end
lower=interpolatedQuantile(x,0.025); upper=interpolatedQuantile(x,0.975);
end

function value=interpolatedQuantile(x,p)
position=1+(numel(x)-1)*p; lo=floor(position); hi=ceil(position);
if lo==hi, value=x(lo); else, value=x(lo)+(position-lo)*(x(hi)-x(lo)); end
end

function [b,V,info]=fitBinomial(X,y,cluster)
y=double(y(:)); assert(all(isfinite(y)&y>=0&y<=1));
[n,nParameters]=size(X); b=zeros(nParameters,1); ridge=1e-8;
converged=false; iteration=0;
for iteration=1:100
  eta=X*b; mu=logistic(eta); w=max(mu.*(1-mu),1e-6);
  Xw=spdiags(sqrt(w),0,n,n)*X;
  H=full(Xw'*Xw)+ridge.*eye(nParameters);
  delta=H\full(X'*(y-mu)); oldNLL=binomialNLL(y,mu); scale=1;
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
info=struct('converged',converged,'iterations',iteration,'nTrials',n, ...
  'nParameters',nParameters,'nClusters',numel(groups), ...
  'logLikelihood',-binomialNLL(y,mu));
end

function p=logistic(x)
p=1./(1+exp(-max(min(x,35),-35)));
end

function value=binomialNLL(y,p)
p=min(max(p,1e-12),1-1e-12);
value=-sum(y.*log(p)+(1-y).*log(1-p));
end
