function tests=testCorrelationLengthProfile
tests=functiontests(localfunctions);
end

function testLikelihoodOnlyContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
idq=string(cache.trials.domain)=="IDQ";
cache.trials.driftDirectionDeg(idq)=30;
sessions=unique(string(cache.trials.sessionID(idq)),'stable');
q=find(string(cache.trials.sessionID)==sessions(1));
cache.trials.driftDirectionDeg(q)=[10;130];
[~,trialGroup]=maxpoolfit.driftRegimeMetadata(cache,'ReadHeaders',false);
adaptationByGroup=table( ...
  ["IDQ-fixed";"IDQ-randomized";"IDS";"IDR"], ...
  [0.1;0.4;0.4;0.6], ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
baseFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
P=maxpoolfit.profileNoiseCorrelationLength(cache,baseFit, ...
  'CorrelationLengthGridDeg',[15 30],'NDraws',8, ...
  'ComputeGainDiagnostics',false);
verifyEqual(testCase,height(P.overallProfile),2);
verifyEqual(testCase,height(P.groupProfile),8);
verifyEqual(testCase,P.correlationLengthGridDeg,[15;30]);
verifyEqual(testCase,max(P.overallProfile.deltaLogLFromBest),0, ...
  'AbsTol',1e-12);
verifyTrue(testCase,all(isfinite(P.groupProfile.noiseSDPC)));
verifyTrue(testCase,all(isnan(P.groupProfile.probeGainRMSE)));
end
