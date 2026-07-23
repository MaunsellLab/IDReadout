function tests=testMTDOGCore
tests=functiontests(localfunctions);
end

function testZeroSurroundIsExactGaussian(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
A=makeGaussianReadoutBank([0 90],5,mt);
B=makeGaussianReadoutBank([0 90],5,mt, ...
  'SurroundMagnitude',0,'SigmaSurroundDeg',nan);
verifyEqual(testCase,A.weightsPhi,B.weightsPhi,'AbsTol',0);
end

function testFiniteSurroundSurvivesMeanSubtraction(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
directions=(0:2:358)';inputs=(0:10:350)';
G=makeGaussianReadoutBank(directions,5,mt);
D=makeGaussianReadoutBank(directions,5,mt, ...
  'SurroundMagnitude',.3,'SigmaSurroundDeg',40);
T=mtPopulationTemplate(inputs,mt);KG=T*G.weightsPhi';KD=T*D.weightsPhi';
scale=sum(KD(:).*KG(:))/sum(KG(:).^2);
verifyGreaterThan(testCase,max(abs(KD-scale.*KG),[],'all'),1e-3);
end

function testDOGRunIDFormatting(testCase)
verifyEqual(testCase,formatDOGRunID("DOG001",5,nan,0), ...
  "DOG001_sigmaC5_baseline");
verifyEqual(testCase,formatDOGRunID("DOG001",5,40,.15), ...
  "DOG001_sigmaC5_sigmaS40_a0p15");
end

function testDOGRunIDRejectsNarrowSurround(testCase)
verifyError(testCase,@()formatDOGRunID("DOG001",5,5,.2), ...
  'formatDOGRunID:BadSurround');
end
