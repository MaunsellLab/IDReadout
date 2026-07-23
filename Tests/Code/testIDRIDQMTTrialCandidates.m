function tests=testIDRIDQMTTrialCandidates
tests=functiontests(localfunctions);
end

function testOpponentSignsAndEqualPeakReadouts(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
condition=toyCondition();
R=evaluateMTTrialCandidateCondition(condition,mt,20,false);
i0=find(R.candidateDirectionsDeg==0); i180=find(R.candidateDirectionsDeg==180);
verifyGreaterThan(testCase,R.change.evidence(1,i0),0);
verifyLessThan(testCase,R.change.evidence(1,i180),0);
verifyLessThan(testCase,R.change.evidence(2,i0),0);
verifyGreaterThan(testCase,R.change.evidence(2,i180),0);
B=makeGaussianReadoutBank([0 180],20,mt);
verifyEqual(testCase,max(B.weightsPhi,[],2),[1;1],'AbsTol',1e-12);
end

function testDuplicateOpponentMechanismsCollapse(testCase)
mt=makeMTReadoutForwardModel();
C=toyCondition(); C.componentDirectionsDeg=[0 180];
C.componentRoles=["primary","probe"];
C.changeComponentCoherence=[1 0;-1 0];
C.noChangeComponentCoherence=zeros(2,2);
C.changeSignalCoherence=zeros(2,2);C.noChangeSignalCoherence=zeros(2,2);
R=evaluateMTTrialCandidateCondition(C,mt,20,true);
verifyEqual(testCase,numel(R.candidateDirectionsDeg),2);
verifyTrue(testCase,all(R.isPhysicalCandidate));
verifyFalse(testCase,any(R.isOpponentOnlyCandidate));
verifyTrue(testCase,contains(R.roleLabels(R.candidateDirectionsDeg==0),'opponentOf_probe'));
end

function testRectifiedActivationsAreNonnegative(testCase)
mt=makeMTReadoutForwardModel();
R=evaluateMTTrialCandidateCondition(toyCondition(),mt,20,true);
verifyGreaterThanOrEqual(testCase,min(R.change.activation,[],'all'),0);
verifyGreaterThanOrEqual(testCase,min(R.noChange.activation,[],'all'),0);
end

function C=toyCondition()
C=struct();C.dataset="toy";C.offsetDeg=120;
C.componentDirectionsDeg=[0 120 -120];
C.componentRoles=["primary","probePlus","probeMinus"];
C.sessionIndex=[1;1];C.trialIndex=[1;2];
C.nExcludedNonfinite=0;
C.pedestalPC=[];
C.changeComponentCoherence=[1 0 0;-1 0 0];
C.noChangeComponentCoherence=zeros(2,3);
C.changeSignalCoherence=zeros(2,3);C.noChangeSignalCoherence=zeros(2,3);
end
