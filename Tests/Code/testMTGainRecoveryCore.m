function tests=testMTGainRecoveryCore
tests=functiontests(localfunctions);
end

function testMatchesAndRestoresSyntheticOrder(testCase)
row=matchSyntheticTrialRows([1;1;2],[10;11;3],[],[2;1],[3;10],[]);
verifyEqual(testCase,row,[3;1]);
end

function testOffsetDisambiguatesTrials(testCase)
row=matchSyntheticTrialRows([1;1],[7;7],[10;25],[1;1],[7;7],[25;10]);
verifyEqual(testCase,row,[2;1]);
end

function testDuplicateSourceRejected(testCase)
verifyError(testCase,@() matchSyntheticTrialRows( ...
  [1;1],[7;7],[],1,7,[]),'matchSyntheticTrialRows:DuplicateSourceKey');
end

function testMissingTrialRejected(testCase)
verifyError(testCase,@() matchSyntheticTrialRows( ...
  1,7,[],1,8,[]),'matchSyntheticTrialRows:MissingTrial');
end
