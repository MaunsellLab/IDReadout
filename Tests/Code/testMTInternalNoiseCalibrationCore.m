function tests=testMTInternalNoiseCalibrationCore
tests=functiontests(localfunctions);
end

function testDenseBankTilesCircle(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:2:358)',20,mt);
verifyEqual(testCase,size(B.weightsPhi,1),180);
verifyEqual(testCase,max(B.weightsPhi,[],2),ones(180,1),'AbsTol',1e-12);
end

function testProjectedCovariance(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:10:350)',20,mt);
N=factorGaussianReadoutNoise(B,'RelativeEigenvalueTolerance',1e-12);
reconstructed=N.factor*N.factor';
verifyEqual(testCase,reconstructed,N.covariance,'RelTol',1e-9,'AbsTol',1e-9);
verifyGreaterThan(testCase,N.covariance(1,2),N.covariance(1,19));
end

function testDenseEvidenceMatchesPopulationCalculation(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:5:355)',20,mt);
coh=[7 7 7;26 7 7];dirs=[0 120 -120];
R=evaluateDenseMTStepCandidates(dirs,coh,B,mt,'RectifyCandidates',false);
M=evaluateMTCoherenceDistribution(dirs,coh,mt);
verifyEqual(testCase,R.evidence,M*B.weightsPhi','AbsTol',1e-10);
end

function testNearbyReadoutsShareMoreNoise(testCase)
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:2:358)',20,mt);
N=factorGaussianReadoutNoise(B);
i0=find(B.candidateDirectionsDeg==0);i2=find(B.candidateDirectionsDeg==2);
i180=find(abs(abs(B.candidateDirectionsDeg)-180)<1e-12,1,'first');
verifyNotEmpty(testCase,i180);
corr02=N.covariance(i0,i2)/sqrt(N.covariance(i0,i0)*N.covariance(i2,i2));
corr180=N.covariance(i0,i180)/sqrt(N.covariance(i0,i0)*N.covariance(i180,i180));
verifyGreaterThan(testCase,corr02,.99);
verifyLessThan(testCase,corr180,.01);
end
