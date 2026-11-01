function tests=testJointGainFit
tests=functiontests(localfunctions);
end

function testJointGainOptimizerContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',20);
sessions=unique(string(cache.trials.sessionID),'stable');
correct=false(height(cache.trials),1);
for k=1:numel(sessions)
  q=string(cache.trials.sessionID)==sessions(k);
  correct(q)=cache.trials.stepMagnitudePC(q)>= ...
    median(cache.trials.stepMagnitudePC(q));
end
cache.trials.correct=correct; cache.trials.choiceRDK=cache.trials.changeRDK;
cache.trials.choiceRDK(~correct)=uint8(3)-cache.trials.changeRDK(~correct);
trialGroup=string(cache.trials.domain); groups=unique(trialGroup,'stable');
adaptationByGroup=table(groups,repmat(.4,numel(groups),1), ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=64;
baseFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
P=maxpoolfit.optimizeGainMechanism(cache,baseFit, ...
  'TargetGroups',["IDQ" "IDS"], ...
  'DirectionalWidthBoundsDeg',[15 30],'AdaptationBounds',[.1 .7], ...
  'CorrelationLengthBoundsDeg',[20 100],'InitialPoint',[20 .4 60], ...
  'NInitialPoints',2,'NStarts',1,'MaxFunctionEvaluationsPerStart',3, ...
  'MinimumGainSessions',1,'NDraws',8,'FinalNDraws',8);
verifyGreaterThanOrEqual(testCase,height(P.evaluation),2);
verifyTrue(testCase,all(isfinite(P.evaluation.changedGainRMSE)));
verifyTrue(testCase,any(P.evaluation.isPareto));
verifyEqual(testCase,P.gainOptimalFit.config.sigmaMTDeg, ...
  P.gainOptimalFit.config.sigmaAdaptDeg);
verifyEqual(testCase,sort(P.gainOptimalFit.groups),sort(["IDQ";"IDS"]));
verifyTrue(testCase,all(P.gainOptimalPoint>=P.parameterLowerBounds));
verifyTrue(testCase,all(P.gainOptimalPoint<=P.parameterUpperBounds));
end
