function tests=testDirectionalWidthGainFit
tests=functiontests(localfunctions);
end

function testGainTargetedProfileContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',20);
% The generic synthetic cache does not impose a psychometric relationship.
% Construct one within every session so the empirical gain normalization has
% the positive reference slope required of real experimental data.
sessions=unique(string(cache.trials.sessionID),'stable');
correct=false(height(cache.trials),1);
for k=1:numel(sessions)
  q=string(cache.trials.sessionID)==sessions(k);
  correct(q)=cache.trials.stepMagnitudePC(q)>= ...
    median(cache.trials.stepMagnitudePC(q));
end
cache.trials.correct=correct;
cache.trials.choiceRDK=cache.trials.changeRDK;
cache.trials.choiceRDK(~correct)=uint8(3)-cache.trials.changeRDK(~correct);
trialGroup=string(cache.trials.domain);
groups=unique(trialGroup,'stable');
adaptationByGroup=table(groups,repmat(.4,numel(groups),1), ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=64;
baseFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
P=maxpoolfit.profileDirectionalWidthByGain(cache,baseFit, ...
  'DirectionalWidthGridDeg',[25 37.5], ...
  'TargetGroups',["IDQ" "IDS"],'MinimumGainSessions',1,'NDraws',8);
verifyEqual(testCase,height(P.overallProfile),2);
verifyEqual(testCase,height(P.groupProfile),6);
verifyEqual(testCase,P.directionalWidthGridDeg,[25;37.5]);
verifyTrue(testCase,all(isfinite(P.overallProfile.changedGainRMSE)));
verifyTrue(testCase,ismember(P.gainOptimalWidthDeg,[25;37.5]));
verifyTrue(testCase,ismember(P.choiceOptimalWidthDeg,[25;37.5]));
verifyEqual(testCase,P.gainOptimalFit.config.sigmaMTDeg, ...
  P.gainOptimalFit.config.sigmaAdaptDeg);
verifyTrue(testCase,all(ismember(unique(string(P.gainProfile.source)), ...
  ["observed" "directional-width model"])));
end
