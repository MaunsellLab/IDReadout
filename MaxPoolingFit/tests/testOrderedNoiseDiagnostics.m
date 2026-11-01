function tests=testOrderedNoiseDiagnostics
tests=functiontests(localfunctions);
end

function testDirectionalChangeGainContract(testCase)
[cache,mask,cfg]=exactConditionCache();
fit=maxpoolfit.fitModel(cache,cfg,'TrialMask',mask,'NDraws',8);
prediction=maxpoolfit.predict(cache,fit,cfg,'TrialMask',mask,'NDraws',8);
result=struct('finalFit',fit,'prediction',prediction);
D=maxpoolfit.decisionRuleGainDiagnostic(cache,result, ...
  'Domain',"IDR",'Animal',"Meetz",'ProbeOffsetDeg',25, ...
  'StepMagnitudePC',15,'FitDraws',8,'PredictionDraws',8, ...
  'MinimumTrials',20,'MinimumSessions',2);
verifyEqual(testCase,D.gainSummary.source, ...
  ["observed";"absolute-step-max";"directional-change-max"]);
verifyEqual(testCase,D.modelComparison.model, ...
  ["absolute-step-max";"directional-change-max"]);
verifyEqual(testCase,D.modelComparison.adaptation(2),0);
verifyEqual(testCase,D.directionalFit.config.decisionRule, ...
  "directional-change-max");
verifyEqual(testCase,height(D.domainComparison),2);
end

function testSharedRDKNoiseSweepContract(testCase)
[cache,mask,cfg]=exactConditionCache();
fit=maxpoolfit.fitModel(cache,cfg,'TrialMask',mask,'NDraws',8);
D=maxpoolfit.sharedRDKNoiseSweep(cache,fit, ...
  'Domain',"IDR",'Animal',"Meetz",'ProbeOffsetDeg',25, ...
  'StepMagnitudePC',15,'NoiseSDGridPC',[0 10], ...
  'RDKNoiseCorrelationGrid',[0 1],'NDraws',8, ...
  'MinimumTrials',20,'MinimumSessions',2);
verifyEqual(testCase,height(D.sweep),4);
verifyEqual(testCase,unique(D.sweep.rdkNoiseCorrelation,'stable'),[0;1]);
verifyEqual(testCase,D.selection.rdkNoiseCorrelationGrid,[0;1]);
verifyTrue(testCase,D.pairedDrawContract);
verifyEqual(testCase,height(D.targets),2);
verifyTrue(testCase,all(isfinite(D.sweep.upstreamAccuracy)));
verifyTrue(testCase,any(D.sweep.downstreamFeasible));
end

function [cache,mask,cfg]=exactConditionCache()
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(863,'twister'); cfg=maxpoolfit.defaultConfig();
cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8; cfg.batchSize=16; cfg.verbose=false;
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
  cache.trials.sessionID(t)="ordered_diagnostic_session_"+(mod(k-1,nSessions)+1);
  if rand<0.75
    cache.trials.choiceRDK(t)=cache.trials.changeRDK(t);
  else
    cache.trials.choiceRDK(t)=uint8(3-double(cache.trials.changeRDK(t)));
  end
end
end
