function tests=testMTGaussianForwardCore
tests=functiontests(localfunctions);
end

function testTemplatesAreMeanZero(testCase)
mt=makeMTReadoutForwardModel();
T=mtPopulationTemplate([0 37 180],mt);
verifyLessThan(testCase,max(abs(mean(T,2))),1e-12);
end

function testDistributionIsLinear(testCase)
mt=makeMTReadoutForwardModel();
M=evaluateMTCoherenceDistribution([0 120],[2 -3],mt);
expected=2*mtPopulationTemplate(0,mt)-3*mtPopulationTemplate(120,mt);
verifyEqual(testCase,M,expected,'AbsTol',1e-12);
end

function testSingleNeuronPreferredNullRatio(testCase)
mt=makeMTReadoutForwardModel();
idx=find(mt.phiDeg==0,1,'first');
preferred=mtPopulationTemplate(0,mt);
null=mtPopulationTemplate(180,mt);
ratio=preferred(idx)/abs(null(idx));
verifyGreaterThan(testCase,ratio,2);
verifyLessThan(testCase,ratio,3);
end

function testGaussianBankRotation(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank([0 90],20,mt);
verifyEqual(testCase,max(B.weightsPhi,[],2),ones(2,1),'AbsTol',1e-12);
verifyGreaterThanOrEqual(testCase,min(B.weightsPhi,[],'all'),0);
end

function testRectificationOccursAfterReadout(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank([0 180],20,mt);
M=evaluateMTCoherenceDistribution(0,-1,mt);
signed=evaluateGaussianDirectionCandidates(M,B,'PoolingMode','max');
rectified=evaluateGaussianDirectionCandidates(M,B, ...
  'RectifyCandidates',true,'PoolingMode','max');
verifyEqual(testCase,rectified.activation,max(signed.evidence,0),'AbsTol',1e-12);
verifyGreaterThanOrEqual(testCase,rectified.pooled,0);
end

function testPNormRequiresNonnegativeCandidates(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank([0 180],20,mt);
M=evaluateMTCoherenceDistribution(0,1,mt);
verifyError(testCase,@() evaluateGaussianDirectionCandidates( ...
  M,B,'PoolingMode','pnorm'),'evaluateGaussianDirectionCandidates:SignedPNorm');
end

function testNegativeAsymptoteIsNestedAndUnitPeak(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
B0=makeGaussianReadoutBank(0,5,mt);
B=makeGaussianReadoutBank(0,5,mt,'NegativeAsymptoteMagnitude',0.2);
verifyEqual(testCase,max(B.weightsPhi),1,'AbsTol',1e-12);
verifyEqual(testCase,B.weightsPhi,1.2.*B0.weightsPhi-0.2,'AbsTol',1e-12);
verifyEqual(testCase,B.negativeAsymptoteMagnitude,0.2);
[~,i180]=min(abs(abs(double(mt.phiDeg))-180));
verifyEqual(testCase,B.weightsPhi(i180),-0.2,'AbsTol',1e-10);
end

function testNegativeAsymptoteRejectsNegativeMagnitude(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
verifyError(testCase,@()makeGaussianReadoutBank(0,5,mt, ...
  'NegativeAsymptoteMagnitude',-0.1), ...
  'MATLAB:InputParser:ArgumentFailedValidation');
end

function testDOGIsNestedAndUnitPeak(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
G=makeGaussianReadoutBank(0,5,mt);
B=makeGaussianReadoutBank(0,5,mt, ...
  'SurroundMagnitude',0.3,'SigmaSurroundDeg',40);
phi=double(mt.phiDeg(:)');distance=mod(phi+180,360)-180;
surround=exp(-(distance.^2)./(2*40^2));
verifyEqual(testCase,B.weightsPhi,1.3.*G.weightsPhi-0.3.*surround, ...
  'AbsTol',1e-12);
verifyEqual(testCase,max(B.weightsPhi),1,'AbsTol',1e-12);
verifyLessThan(testCase,min(B.weightsPhi),0);
end

function testDOGRequiresBroaderSurround(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
verifyError(testCase,@()makeGaussianReadoutBank(0,5,mt, ...
  'SurroundMagnitude',0.2,'SigmaSurroundDeg',5), ...
  'makeGaussianReadoutBank:BadSurroundSigma');
end
