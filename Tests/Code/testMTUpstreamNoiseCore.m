function tests=testMTUpstreamNoiseCore
tests=functiontests(localfunctions);
end

function testCandidateCovarianceFactor(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:10:350)',20,mt);
N=factorUpstreamDirectionalNoise((0:359)',B,mt, ...
  'RelativeEigenvalueTolerance',1e-12);
verifyEqual(testCase,N.factor*N.factor',N.candidateCovariance, ...
  'RelTol',1e-9,'AbsTol',1e-8);
end

function testTuningCreatesExpectedCorrelationSigns(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
B=makeGaussianReadoutBank((0:10:350)',20,mt);
N=factorUpstreamDirectionalNoise((0:359)',B,mt);
phi=mt.phiDeg(:);i0=find(phi==0);i10=find(phi==10);
i60=find(phi==60);i90=find(phi==90);
i180=find(abs(abs(phi)-180)<1e-12,1,'first');
verifyGreaterThan(testCase,N.mtCorrelation(i0,i10),.9);
verifyGreaterThan(testCase,N.mtCorrelation(i0,i60),0);
verifyLessThan(testCase,N.mtCorrelation(i0,i90),0);
verifyLessThan(testCase,N.mtCorrelation(i0,i180),-.4);
end

function testNoSeparateReadoutNoise(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:5:355)',20,mt);
N=factorUpstreamDirectionalNoise((0:359)',B,mt);
u=zeros(1,360);u(1)=1;
candidateFromFactor=u*N.inputToCandidate;
mtResponse=u*N.inputToMT;
candidateFromMT=mtResponse*B.weightsPhi';
verifyEqual(testCase,candidateFromFactor,candidateFromMT,'AbsTol',1e-12);
end

function testTrialInputsRetainIDQPreStepPedestalField(testCase)
% Structural guard: the production input builder must retain the new epoch.
source=fileread(which('makeIDRIDQMTTrialInputs'));
verifyTrue(testCase,contains(source,'preStepComponentCoherence'));
end
