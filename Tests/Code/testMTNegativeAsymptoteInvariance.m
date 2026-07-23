function tests=testMTNegativeAsymptoteInvariance
tests=functiontests(localfunctions);
end

function testConstantTailOnlyScalesCandidateGeometry(testCase)
S=analyzeMTNegativeAsymptoteInvariance( ...
  'NegativeAsymptoteMagnitude',[0 .05 .2 .5]);
verifyLessThan(testCase,max(S.maxAbsTemplateSum),1e-10);
verifyLessThan(testCase,max(S.maxAbsSignalScaleResidual),1e-9);
verifyLessThan(testCase,max(S.maxAbsCovarianceScaleResidual),1e-6);
verifyEqual(testCase,S.choiceAgreement,ones(height(S),1));
end

function testZeroAsymptoteIsOriginalBank(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
A=makeGaussianReadoutBank([0 90],5,mt);
B=makeGaussianReadoutBank([0 90],5,mt, ...
  'NegativeAsymptoteMagnitude',0);
verifyEqual(testCase,A.weightsPhi,B.weightsPhi,'AbsTol',0);
end
