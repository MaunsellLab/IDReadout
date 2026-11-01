function tests=testOpponentAxisDecomposition
tests=functiontests(localfunctions);
end

function testOpponentAxisContract(testCase)
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
fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup,adaptationByGroup, ...
  'NDraws',8);
joint=struct('gainOptimalFit',fit,'choiceOptimalFit',fit);
P=maxpoolfit.opponentAxisDecomposition(cache,joint, ...
  'ProbeOffsetDeg',45,'FitSelections',"gain-optimal",'NDraws',8);
verifyEqual(testCase,height(P.summary),3);
verifyEqual(testCase,P.summary.configuration, ...
  ["paired actual";"single 180 A";"single 180 B"]);
verifyTrue(testCase,all(isfinite(P.summary.meanDriftMinusOpponentResponse)));
verifyTrue(testCase,all(P.summary.internalProbeWinnerFraction>=0 & ...
  P.summary.internalProbeWinnerFraction<=1));
verifyTrue(testCase,all(isfinite(P.summary.predictedChangedProbeGain)));
end
