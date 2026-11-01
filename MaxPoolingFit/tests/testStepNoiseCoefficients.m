function tests=testStepNoiseCoefficients
tests=functiontests(localfunctions);
end

function testExactStepCoefficientContract(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8;
cfg.finalNoiseDraws=8; cfg.batchSize=16; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',40);
% Exact-step diagnostics require repeated step levels. syntheticCache uses a
% continuous random step by default, so put every trial at the same exact
% level while preserving its original drift-noise residual.
fixedStep=15;
for k=1:height(cache.trials)
  delta=fixedStep-cache.trials.stepMagnitudePC(k);
  cache.trials.stepMagnitudePC(k)=fixedStep;
  cache.trials.stepEndpointCohPC(k)= ...
    cache.trials.stepEndpointCohPC(k)+delta;
  changed=double(cache.trials.changeRDK(k));
  cache.stepCohPC{k,changed}(1)=cache.stepCohPC{k,changed}(1)+delta;
end
mask=cfg.primarySelector(cache.trials);
p1=nan(height(cache.trials),1); noise=[13 8 10];
domains=unique(cache.trials.domain(mask),'stable');
for k=1:numel(domains)
  idx=find(mask & cache.trials.domain==domains(k));
  p1(idx)=maxpoolfit.choiceProbability(cache,idx,cfg.adaptation,noise(k),cfg, ...
    'NDraws',8);
end
prediction=struct('pRDK1',p1,'trialIdx',find(mask));
D=maxpoolfit.noiseCoefficientsByStep(cache,prediction, ...
  'TrialMask',mask,'MinimumTrials',20);
S=D.summary;
verifyFalse(testCase,isempty(S));
verifyEqual(testCase,unique(S.source,'stable'), ...
  ["observed";"primary domain noise"]);
verifyTrue(testCase,all(ismember( ...
  {'domain','animal','probeOffsetDeg','stepMagnitudePC','nTrials','nSessions', ...
   'accuracy','driftGainImpliedWinFraction','driftGainImpliedWinFractionSE', ...
   'driftGainSignsConsistent','probeGainImpliedWinFraction', ...
   'probeGainImpliedWinFractionSE','probeGainSignsConsistent', ...
   'changedDriftCoefficientPerPC','unchangedDriftCoefficientPerPC', ...
   'changedProbeCoefficientPerPC','unchangedProbeCoefficientPerPC'}, ...
  S.Properties.VariableNames)));
verifyTrue(testCase,all(isfinite(S.changedDriftCoefficientPerPC)));
verifyTrue(testCase,all(isfinite(S.unchangedDriftCoefficientPerPC)));
verifyTrue(testCase,all(isfinite(S.changedProbeCoefficientPerPC)));
verifyTrue(testCase,all(isfinite(S.unchangedProbeCoefficientPerPC)));
verifyTrue(testCase,all(S.accuracy>=0 & S.accuracy<=1));
q=string(S.source)=="primary domain noise" & S.driftGainSignsConsistent;
verifyTrue(testCase,any(q));
expected=S.changedDriftCoefficientPerPC(q)./( ...
  S.changedDriftCoefficientPerPC(q)-S.unchangedDriftCoefficientPerPC(q));
verifyEqual(testCase,S.driftGainImpliedWinFraction(q),expected,'AbsTol',1e-12);
verifyEqual(testCase,D.nTrials,sum(mask));
end

function testSparseExactStepsReturnDefinedEmptyTable(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
mask=cache.trials.valid & cache.trials.hasStepNoise & cache.trials.isIncrement;
p1=nan(height(cache.trials),1); p1(mask)=0.5;
prediction=struct('pRDK1',p1,'trialIdx',find(mask));
D=maxpoolfit.noiseCoefficientsByStep(cache,prediction, ...
  'TrialMask',mask,'MinimumTrials',20);
verifyEqual(testCase,height(D.summary),0);
verifyTrue(testCase,all(ismember({'source','accuracy', ...
  'driftGainImpliedWinFraction','probeGainImpliedWinFraction'}, ...
  D.summary.Properties.VariableNames)));
end
