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

function testRectifiedMTTemplateMatchesHistoricalDefinition(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
S=mtPopulationTemplate([0 45 180],mt,'MTTemplateMode','signed');
R=mtPopulationTemplate([0 45 180],mt,'MTTemplateMode','rectified');
verifyEqual(testCase,R,max(S,0),'AbsTol',0);
verifyGreaterThanOrEqual(testCase,min(R,[],'all'),0);
verifyGreaterThan(testCase,min(sum(R,2)),0);
end

function testPhysicalCoherenceSignSurvivesTemplateRectification(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
B=makeGaussianReadoutBank((0:10:350)',5,mt, ...
  'SurroundMagnitude',.15,'SigmaSurroundDeg',80);
P=evaluateDenseMTStepCandidates(0,1,B,mt,'RectifyCandidates',false, ...
  'MTTemplateMode','rectified');
N=evaluateDenseMTStepCandidates(0,-1,B,mt,'RectifyCandidates',false, ...
  'MTTemplateMode','rectified');
verifyEqual(testCase,N.evidence,-P.evidence,'AbsTol',1e-12);
end

function testRectifiedUpstreamNoiseUsesSameTemplates(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);
inputs=(0:10:350)';B=makeGaussianReadoutBank((0:10:350)',5,mt);
N=factorUpstreamDirectionalNoise(inputs,B,mt, ...
  'MTTemplateMode','rectified');
expected=mtPopulationTemplate(inputs,mt,'MTTemplateMode','rectified')* ...
  B.weightsPhi';
verifyEqual(testCase,N.inputToCandidate,expected,'AbsTol',1e-12);
verifyEqual(testCase,N.mtTemplateMode,"rectified");
end

function testRectificationMakesConstantTailNondegenerate(testCase)
mt=makeMTReadoutForwardModel('sigmaMTDeg',37.5);inputs=(0:10:350)';
G=makeGaussianReadoutBank((0:10:350)',5,mt);
A=makeGaussianReadoutBank((0:10:350)',5,mt, ...
  'NegativeAsymptoteMagnitude',.2);
T=mtPopulationTemplate(inputs,mt,'MTTemplateMode','rectified');
KG=T*G.weightsPhi';KA=T*A.weightsPhi';
scale=sum(KA(:).*KG(:))/sum(KG(:).^2);
verifyGreaterThan(testCase,max(abs(KA-scale.*KG),[],'all'),1e-3);
end
