function tests=testMTReadoutWidthSweepCore
tests=functiontests(localfunctions);
end

function testIntegerWidthToken(testCase)
verifyEqual(testCase,formatReadoutWidthRunID("Sweep001",5), ...
  "Sweep001_sigmaR5");
end

function testFractionalWidthToken(testCase)
verifyEqual(testCase,formatReadoutWidthRunID("Sweep001",2.5), ...
  "Sweep001_sigmaR2p5");
end

function testBadWidthRejected(testCase)
verifyError(testCase,@() formatReadoutWidthRunID("Sweep001",0), ...
  'formatReadoutWidthRunID:BadWidth');
end

function testLongPrefixRejected(testCase)
verifyError(testCase,@() formatReadoutWidthRunID(repmat('a',1,78),20), ...
  'formatReadoutWidthRunID:RunIDTooLong');
end
