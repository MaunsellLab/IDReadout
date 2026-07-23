function tests=testMTUpstreamChoiceSimulation
tests=functiontests(localfunctions);
end

function testZeroNoiseSelectsChangedPatch(testCase)
[C,B,N,mt]=toyInputs();stream=RandStream('mt19937ar','Seed',11);
R=simulateUpstreamChoiceCondition(C,B,N,mt,zeros(4,1),3,stream);
verifyTrue(testCase,all(R.correctReplicates,'all'));
verifyEqual(testCase,R.trialProbabilityCorrect,ones(4,1));
end

function testReproducibleChoices(testCase)
[C,B,N,mt]=toyInputs();sigma=ones(4,1)*2;
s1=RandStream('mt19937ar','Seed',17);s2=RandStream('mt19937ar','Seed',17);
R1=simulateUpstreamChoiceCondition(C,B,N,mt,sigma,5,s1);
R2=simulateUpstreamChoiceCondition(C,B,N,mt,sigma,5,s2);
verifyEqual(testCase,R1.correctReplicates,R2.correctReplicates);
verifyEqual(testCase,R1.changeWinnerFraction,R2.changeWinnerFraction);
end

function testChoiceStorageIsCompact(testCase)
[C,B,N,mt]=toyInputs();stream=RandStream('mt19937ar','Seed',23);
R=simulateUpstreamChoiceCondition(C,B,N,mt,ones(4,1),7,stream);
verifyClass(testCase,R.correctReplicates,'logical');
verifyFalse(testCase,isfield(R,'candidateEvidence'));
end

function [C,B,N,mt]=toyInputs()
mt=makeMTReadoutForwardModel();
B=makeGaussianReadoutBank((0:10:350)',20,mt);
N=factorUpstreamDirectionalNoise((0:10:350)',B,mt);
C=struct('dataset',"toy",'offsetDeg',120, ...
  'componentDirectionsDeg',[0 120 -120], ...
  'changeComponentCoherence',repmat([20 0 0],4,1), ...
  'noChangeComponentCoherence',zeros(4,3), ...
  'sessionIndex',ones(4,1),'trialIndex',(1:4)', ...
  'experimentalCorrect',true(4,1));
end
