function tests=testPerformanceMatchedFit
tests=functiontests(localfunctions);
end

function testPerformanceMatchedOutputContract(testCase)
cfg=maxpoolfit.defaultConfig();
cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[fit,prediction]=maxpoolfit.fitPerformanceMatched(cache,cfg, ...
  'NDraws',8,'GroupVariables',{'domain'},'CalibrationGridSize',5);
verifyEqual(testCase,fit.adaptation,cfg.adaptation);
verifyEqual(testCase,fit.groupVariables,"domain");
verifyEqual(testCase,height(fit.groupFit),3);
verifyEqual(testCase,sum(fit.groupFit.nTrials),fit.nTrials);
verifyTrue(testCase,all(isfinite(fit.groupFit.noiseSDPC)));
verifyTrue(testCase,all(isfinite(prediction.pRDK1(prediction.trialIdx))));
verifyEqual(testCase,prediction.trialIdx,find(fit.trialMask));
end

function testPerformanceMatchedAcceptsGroupAdaptation(testCase)
cfg=maxpoolfit.defaultConfig();
cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
trialAdaptation=nan(height(cache.trials),1);
trialAdaptation(string(cache.trials.domain)=="IDR")=0.6;
trialAdaptation(string(cache.trials.domain)=="IDQ")=0.4;
trialAdaptation(string(cache.trials.domain)=="IDS")=0.4;
[fit,prediction]=maxpoolfit.fitPerformanceMatched(cache,cfg, ...
  'AdaptationByTrial',trialAdaptation,'NDraws',8, ...
  'GroupVariables',{'domain'},'CalibrationGridSize',5);
verifyTrue(testCase,isnan(fit.adaptation));
verifyEqual(testCase,fit.trialAdaptation,trialAdaptation);
verifyEqual(testCase,prediction.trialAdaptation,trialAdaptation);
verifyTrue(testCase,all(isfinite(prediction.pRDK1(prediction.trialIdx))));
end
