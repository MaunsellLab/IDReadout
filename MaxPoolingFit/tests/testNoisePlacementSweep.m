function tests=testNoisePlacementSweep
tests=functiontests(localfunctions);
end

function testNoisePlacementContract(testCase)
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(761,'twister');
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8;
cfg.batchSize=16; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',80);
mask=cache.trials.domain=="IDR" & cache.trials.animal=="Meetz" & ...
  cache.trials.probeOffsetDeg==25;
idx=find(mask); fixedStep=15; nSessions=8;
for k=1:numel(idx)
  t=idx(k); delta=fixedStep-cache.trials.stepMagnitudePC(t);
  cache.trials.stepMagnitudePC(t)=fixedStep;
  cache.trials.stepEndpointCohPC(t)=cache.trials.stepEndpointCohPC(t)+delta;
  changed=double(cache.trials.changeRDK(t));
  cache.stepCohPC{t,changed}(1)=cache.stepCohPC{t,changed}(1)+delta;
  cache.trials.sessionID(t)="noise_placement_session_"+(mod(k-1,nSessions)+1);
  if rand<0.75
    cache.trials.choiceRDK(t)=cache.trials.changeRDK(t);
  else
    cache.trials.choiceRDK(t)=uint8(3-double(cache.trials.changeRDK(t)));
  end
end
fit=struct('adaptation',cfg.adaptation,'config',cfg,'trialMask',mask, ...
  'modelVersion',cfg.modelVersion,'domains',"IDR",'noiseSDPC',10, ...
  'domainFit',table("IDR",sum(mask),10,0,0, ...
  'VariableNames',{'domain','nTrials','noiseSDPC','logLikelihood','evalCount'}));

D=maxpoolfit.noisePlacementSweep(cache,fit,'Domain',"IDR",'Animal',"Meetz", ...
  'ProbeOffsetDeg',25,'StepMagnitudePC',fixedStep, ...
  'NoiseSDGridPC',[0 5 10],'NDraws',8, ...
  'MinimumTrials',20,'MinimumSessions',2);

verifyEqual(testCase,height(D.sweep),3);
verifyEqual(testCase,D.sweep.noiseSDPC,[0;5;10]);
verifyEqual(testCase,D.selection.stepMagnitudePC,fixedStep);
verifyEqual(testCase,height(D.targets),2);
verifyEqual(testCase,D.targets.component,["drift";"probe"]);
verifyTrue(testCase,all(ismember({'requiredFlipRate','downstreamFeasible', ...
  'upstreamDriftUnchangedOverChanged','postFlipDriftUnchangedOverChanged', ...
  'upstreamProbeUnchangedOverChanged','postFlipProbeUnchangedOverChanged', ...
  'upstreamRatioRMSE','postFlipRatioRMSE','ratioRMSE', ...
  'postFlipAbsoluteGainRMSE'},D.sweep.Properties.VariableNames)));
q=D.sweep.downstreamFeasible;
verifyTrue(testCase,any(q));
verifyEqual(testCase,D.sweep.postFlipAccuracy(q), ...
  repmat(D.observedAccuracy,sum(q),1),'AbsTol',1e-10);
verifyTrue(testCase,all(D.sweep.requiredFlipRate(q)>=0 & ...
  D.sweep.requiredFlipRate(q)<=0.5));
end
