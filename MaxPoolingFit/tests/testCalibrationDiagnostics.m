function tests=testCalibrationDiagnostics
tests=functiontests(localfunctions);
end

function testCalibrationTablesUseSelectedTrials(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
n=height(cache.trials);
mask=false(n,1); mask(1:2:n)=true;
idx=find(mask);
p1=0.25+0.5.*double(cache.trials.changeRDK==1);
p1=min(max(p1,1e-4),1-1e-4);
prediction=struct('trialIdx',idx,'pRDK1',p1);
C=maxpoolfit.calibrationDiagnostics(cache,prediction,'NBins',4);
verifyEqual(testCase,C.nTrials,numel(idx));
verifyEqual(testCase,sum(C.groupSummary.nTrials),numel(idx));
verifyEqual(testCase,sum(C.sessionSummary.nTrials),numel(idx));
verifyEqual(testCase,sum(C.reliability.nTrials),numel(idx));
verifyTrue(testCase,all(isfinite(C.groupSummary.logLoss)));
verifyTrue(testCase,all(isfinite(C.groupSummary.calibrationIntercept)));
verifyTrue(testCase,all(isfinite(C.groupSummary.calibrationSlope)));
end
