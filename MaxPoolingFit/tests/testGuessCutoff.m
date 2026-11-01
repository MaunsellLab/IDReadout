function tests=testGuessCutoff
tests=functiontests(localfunctions);
end

function testNoGuessExactlyMatchesAbsoluteMax(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.batchSize=8;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
idx=1:5; [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
z1=maxpoolfit.noiseDraws(cfg,8,0);
z2=maxpoolfit.noiseDraws(cfg,8,1000003);
p0=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,Rpre,cfg.adaptation,10,cfg,z1,z2,false);
[pGuess,details]=maxpoolfit.choiceProbabilityDetectionGuessFromResponses( ...
  B,A,Rpre,cfg.adaptation,10,-Inf,cfg,z1,z2,0.5);
verifyEqual(testCase,pGuess,p0,'AbsTol',1e-12);
verifyEqual(testCase,details.guessProbability,zeros(numel(idx),1));
end

function testHighCutoffForcesUnbiasedGuessing(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A,Rpre]=maxpoolfit.deterministicResponse(cache,1:5,cfg);
z1=maxpoolfit.noiseDraws(cfg,8,0);
z2=maxpoolfit.noiseDraws(cfg,8,1000003);
[pGuess,details]=maxpoolfit.choiceProbabilityDetectionGuessFromResponses( ...
  B,A,Rpre,cfg.adaptation,10,100,cfg,z1,z2,0.5);
verifyEqual(testCase,pGuess,0.5.*ones(5,1),'AbsTol',1e-12);
verifyEqual(testCase,details.guessProbability,ones(5,1),'AbsTol',1e-12);
end

function testFiniteCutoffOnlyChangesGuessedPairs(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A,Rpre]=maxpoolfit.deterministicResponse(cache,1:5,cfg);
z1=maxpoolfit.noiseDraws(cfg,8,0);
z2=maxpoolfit.noiseDraws(cfg,8,1000003);
p0=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,Rpre,cfg.adaptation,10,cfg,z1,z2,false);
[pCutoff,details]=maxpoolfit.choiceProbabilityDetectionGuessFromResponses( ...
  B,A,Rpre,cfg.adaptation,10,0,cfg,z1,z2,0.5);
verifyLessThanOrEqual(testCase,abs(pCutoff-p0), ...
  details.guessProbability+1e-12);
end

function testGuessCutoffProfileContract(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8;
cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
G=maxpoolfit.profileGuessCutoff(cache,cfg,'NDraws',8, ...
  'KappaGrid',[-Inf 0 1]);
verifyEqual(testCase,G.profile.kappa,[-Inf;0;1]);
verifyEqual(testCase,G.profile.deltaLogLFromNoGuess(1),0,'AbsTol',1e-12);
verifySize(testCase,G.profile.noiseSDPC,[3 3]);
verifySize(testCase,G.profile.guessFraction,[3 3]);
verifyTrue(testCase,all(isfinite(G.profile.logLikelihood)));
verifyTrue(testCase,all(G.profile.guessFraction(1,:)==0));
verifyTrue(testCase,all(G.profile.guessFraction(2:end,:)>=0 & ...
  G.profile.guessFraction(2:end,:)<=1,'all'));
verifyTrue(testCase,istable(G.conditionSummary));
verifyTrue(testCase,all(isfinite(G.conditionSummary.meanGuessProbability)));
end
