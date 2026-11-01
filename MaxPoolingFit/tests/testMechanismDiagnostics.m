function tests=testMechanismDiagnostics
tests=functiontests(localfunctions);
end

function testMechanismDiagnosticOutput(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.finalNoiseDraws=8;
cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
mask=cfg.primarySelector(cache.trials);
primaryFit=struct('adaptation',cfg.adaptation, ...
  'domains',unique(cache.trials.domain(mask),'stable'), ...
  'noiseSDPC',[8;8;8]);
[performanceFit,~]=maxpoolfit.fitPerformanceMatched(cache,cfg, ...
  'TrialMask',mask,'NDraws',8,'GroupVariables',{'domain'}, ...
  'CalibrationGridSize',5);
performanceFit.groupFit.primaryPredictedAccuracy= ...
  performanceFit.groupFit.predictedAccuracy;
D=maxpoolfit.mechanismDiagnostics(cache,primaryFit,performanceFit,1, ...
  'NDraws',8,'MaxTrialsPerCondition',5,'NoiseGridPC',[0 8]);
verifyEqual(testCase,height(D.conditionSummary),1);
verifyGreaterThanOrEqual(testCase,height(D.noiseSweep),2);
verifyEqual(testCase,sum(D.winnerHistogram.changedWinnerFraction),1,'AbsTol',1e-12);
verifyTrue(testCase,isfinite(D.conditionSummary.meanBetweenRDKMargin));
end
