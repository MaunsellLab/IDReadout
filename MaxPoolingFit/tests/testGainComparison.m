function tests=testGainComparison
tests=functiontests(localfunctions);
end

function testPredictorReconstruction(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
idx=find(string(cache.trials.domain)=="IDR",1);
P=maxpoolfit.gainPredictors(cache,idx);
c=double(cache.trials.changeRDK(idx)); u=3-c;
verifyEqual(testCase,P.changedDriftNoisePC, ...
  cache.stepCohPC{idx,c}(1)-cache.trials.stepEndpointCohPC(idx), ...
  'AbsTol',1e-12);
verifyEqual(testCase,P.unchangedDriftNoisePC, ...
  cache.stepCohPC{idx,u}(1)-cache.trials.baseCohPC(idx), ...
  'AbsTol',1e-12);
verifyEqual(testCase,P.changedProbeNoisePC,cache.stepCohPC{idx,c}(2), ...
  'AbsTol',1e-12);
end

function testObservedAndPredictedGainContract(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8;
cfg.batchSize=16; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',80,'NoNoiseFraction',0.5);
noise=[15 10 9];
p1=nan(height(cache.trials),1);
domains=unique(cache.trials.domain,'stable');
for k=1:numel(domains)
  idx=find(cache.trials.domain==domains(k));
  p1(idx)=maxpoolfit.choiceProbability(cache,idx,cfg.adaptation,noise(k),cfg, ...
    'NDraws',8);
end
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>; rng(1441,'twister');
cache.trials.choiceRDK=uint8(1+(rand(height(cache.trials),1)>=p1));
% Make the small synthetic no-noise reference set deterministically
% psychometric. Its purpose is to test gain scaling, not whether a random
% sample happens to produce a positive slope in every domain-animal group.
noNoise=~cache.trials.hasStepNoise;
sessions=unique(cache.trials.sessionID(noNoise),'stable');
for s=1:numel(sessions)
  q=find(noNoise & cache.trials.sessionID==sessions(s));
  isCorrect=cache.trials.stepMagnitudePC(q)>=median(cache.trials.stepMagnitudePC(q));
  chosen=uint8(3-double(cache.trials.changeRDK(q)));
  chosen(isCorrect)=cache.trials.changeRDK(q(isCorrect));
  cache.trials.choiceRDK(q)=chosen;
end
cache.trials.correct=cache.trials.choiceRDK==cache.trials.changeRDK;
mask=cfg.primarySelector(cache.trials);
prediction=struct('pRDK1',p1,'trialIdx',find(mask));
C=maxpoolfit.gainComparison(cache,prediction,prediction,'TrialMask',mask);
verifyEqual(testCase,C.nTrials,sum(mask));
verifyEqual(testCase,unique(C.summary.source,'stable'), ...
  ["observed";"primary domain noise";"performance matched"]);
verifyEqual(testCase,height(C.referencePsychometric),6);
verifyTrue(testCase,all(C.referencePsychometric.slopePerPC>0));
verifyTrue(testCase,all(isfinite(C.summary.changedDriftGain)));
verifyTrue(testCase,all(isfinite(C.summary.changedProbeGain)));
verifyTrue(testCase,ismember('uncertaintyStatus', ...
  C.summary.Properties.VariableNames));

% Exercise the production fallback used when ProbeSessions omit no-noise
% trials: estimate the common scale from noisy trials with offset effects.
fallbackCache=cache; fallbackCache.trials.hasStepNoise(:)=true;
fallbackMask=cfg.primarySelector(fallbackCache.trials);
fallbackPrediction=struct('pRDK1',p1,'trialIdx',find(fallbackMask));
F=maxpoolfit.gainComparison(fallbackCache,fallbackPrediction, ...
  fallbackPrediction,'TrialMask',fallbackMask);
verifyTrue(testCase,all(F.referencePsychometric.method== ...
  "step-noise trials with probe-offset fixed effects"));
verifyTrue(testCase,all(F.referencePsychometric.slopePerPC>0));
end
