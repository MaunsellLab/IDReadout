function tests=testTemporalGainDiagnostics
tests=functiontests(localfunctions);
end

function testFiveBinTemporalContract(testCase)
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(417,'twister');
cache=maxpoolfit.syntheticCache('TrialsPerCondition',100);
n=height(cache.trials); raw=cell(n,2); driver=cell(n,1); p1=nan(n,1);
for t=1:n
  nStreams=numel(cache.streamDirDeg{t,1});
  for j=1:2
    baseline=double(cache.stepCohPC{t,j}(:))'; C=repmat(baseline,30,1);
    for b=1:5
      rows=(b-1)*6+(1:6); noise=2.*randn(1,nStreams);
      if string(cache.trials.domain(t))=="IDR"&&nStreams>=3, noise(3)=noise(2); end
      C(rows,:)=C(rows,:)+noise;
    end
    raw{t,j}=C;
  end
  driver{t}=zeros(30,1);
end
cache.raw=struct('stepStreamCohPC',{raw},'kernelDriverByFrame',{driver});

mask=cache.trials.domain=="IDR" & cache.trials.animal=="Meetz" & ...
  cache.trials.probeOffsetDeg==25;
idx=find(mask); nSessions=8;
for k=1:numel(idx)
  cache.trials.sessionID(idx(k))="temporal_session_"+(mod(k-1,nSessions)+1);
  c=double(cache.trials.changeRDK(idx(k))); u=3-c;
  C=raw{idx(k),c}; U=raw{idx(k),u};
  changed=mean(C(1:6,1))-cache.trials.stepEndpointCohPC(idx(k));
  unchanged=mean(U(1:6,1))-cache.trials.baseCohPC(idx(k));
  pCorrect=1./(1+exp(-(0.7+0.12.*changed-0.06.*unchanged)));
  isCorrect=rand<pCorrect;
  if isCorrect
    cache.trials.choiceRDK(idx(k))=cache.trials.changeRDK(idx(k));
  else
    cache.trials.choiceRDK(idx(k))=uint8(3-double(cache.trials.changeRDK(idx(k))));
  end
  if c==1, p1(idx(k))=pCorrect; else, p1(idx(k))=1-pCorrect; end
end
maxpoolfit.validateCache(cache);
prediction=struct('pRDK1',p1,'trialIdx',idx);
D=maxpoolfit.temporalGainDiagnostics(cache,prediction, ...
  'TrialMask',mask,'Domain',"IDR",'Animal',"Meetz",'ProbeOffsetDeg',25, ...
  'StepMagnitudePC',[],'MinimumTrials',20,'MinimumSessions',2, ...
  'NBootstrap',8,'RandomSeed',19);

verifyEqual(testCase,D.nFrames,30);
verifyEqual(testCase,D.nBins,5);
verifyEqual(testCase,D.binDurationMS,50,'AbsTol',1e-12);
verifyEqual(testCase,height(D.profile),40);
verifyEqual(testCase,height(D.shapeTests),4);
verifyEqual(testCase,unique(D.profile.timeStartMS)',[0 50 100 150 200]);
verifyEqual(testCase,unique(D.profile.timeEndMS)',[50 100 150 200 250]);
verifyEqual(testCase,unique(D.profile.source,'stable'),["observed";"primary domain noise"]);
verifyTrue(testCase,all(ismember({'shapePValue','unchangedOverChanged', ...
  'predictorConditionNumber'},D.shapeTests.Properties.VariableNames)));
verifyTrue(testCase,all(isfinite(D.profile.signedCoefficientPerPC)));
end

function testRejectsNonDivisibleFrameCount(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',20);
n=height(cache.trials); raw=cell(n,2); driver=cell(n,1);
for t=1:n
  for j=1:2, raw{t,j}=repmat(double(cache.stepCohPC{t,j}(:))',29,1); end
  driver{t}=zeros(29,1);
end
cache.raw=struct('stepStreamCohPC',{raw},'kernelDriverByFrame',{driver});
mask=cache.trials.domain=="IDR" & cache.trials.animal=="Meetz" & ...
  cache.trials.probeOffsetDeg==25;
p1=nan(n,1); p1(mask)=0.5; prediction=struct('pRDK1',p1,'trialIdx',find(mask));
verifyError(testCase,@()maxpoolfit.temporalGainDiagnostics( ...
  cache,prediction,'TrialMask',mask,'MinimumTrials',20, ...
  'MinimumSessions',2,'NBootstrap',0), ...
  'maxpoolfit:temporalGainDiagnostics:NonDivisibleFrameCount');
end
