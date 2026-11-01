function tests=testMaxPoolingCore
tests=functiontests(localfunctions);
end

function testCacheAndResponses(testCase)
cfg=maxpoolfit.defaultConfig();
verifyEqual(testCase,cfg.modelVersion, ...
  "DenseHardMax-v5-shared-prestep-adaptation-random-dot-redistribution");
verifyEqual(testCase,cfg.decisionRule,"absolute-step-max");
verifyEqual(testCase,cfg.coherenceResponseMode,"random-dot-redistribution");
verifyEqual(testCase,cfg.adaptation,0.56,'AbsTol',1e-12);
verifyEqual(testCase,cfg.adaptationExponent,0,'AbsTol',1e-12);
verifyEqual(testCase,cfg.adaptationGainStrength,0,'AbsTol',1e-12);
verifyEqual(testCase,cfg.adaptationReferencePC,40,'AbsTol',1e-12);
cache=maxpoolfit.syntheticCache('TrialsPerCondition',2);
maxpoolfit.validateCache(cache);
[B,A,P]=maxpoolfit.deterministicResponse(cache,1:3,cfg);
verifySize(testCase,B,[numel(cfg.directionGridDeg) 2 3]);
verifySize(testCase,A,size(B));
verifySize(testCase,P,size(B));
verifyLessThan(testCase,min(A,[],'all'),0);
verifyLessThan(testCase,min(P,[],'all'),0);
end

function testExponentialAdaptationGainIsNested(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureCoherenceResponse(cfg,"additive");
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A,P]=maxpoolfit.deterministicResponse(cache,1,cfg);
D0=maxpoolfit.decisionResponse(B,A,P,0.56,cfg);
verifyEqual(testCase,D0,B-0.56.*A,'AbsTol',1e-12);

cfg=maxpoolfit.configureAdaptationGain(cfg,0.5,40);
verifyEqual(testCase,cfg.modelVersion, ...
  "DenseHardMax-v14-increment-only-adaptation-gain");
D=maxpoolfit.decisionResponse(B,A,P,0.56,cfg);
expected=P-0.56.*A+exp(-0.5.*P./40).*(B-P);
verifyEqual(testCase,D,expected,'AbsTol',1e-12);
unchanged=maxpoolfit.decisionResponse(P,A,P,0.56,cfg);
verifyEqual(testCase,unchanged,P-0.56.*A,'AbsTol',1e-12);
verifyGreaterThan(testCase,min(maxpoolfit.adaptationGain(P,cfg),[],'all'),0);
end

function testRandomDotRedistributionContract(testCase)
cfg=maxpoolfit.defaultConfig();
theta=cfg.directionGridDeg(:);
mu=maxpoolfit.randomDirectionMean(cfg.sigmaMTDeg);
verifyEqual(testCase,mu,0.261106697653484,'AbsTol',1e-12);
G=maxpoolfit.coherenceTuning(theta,0,cfg.sigmaMTDeg,cfg);
verifyEqual(testCase,mean(G),0,'AbsTol',3e-9);
[~,i45]=min(abs(theta-45));
[~,i120]=min(abs(theta-120));
verifyGreaterThan(testCase,G(i45),0);
verifyLessThan(testCase,G(i120),0);
verifyEqual(testCase,maxpoolfit.adaptationDrive([-1 1],cfg),[-1 1], ...
  'AbsTol',1e-12);
verifyEqual(testCase,maxpoolfit.adaptationGain([-1 1],cfg),[1 1], ...
  'AbsTol',1e-12);

legacy=rmfield(cfg,'coherenceResponseMode');
verifyEqual(testCase,maxpoolfit.coherenceTuning(theta,0, ...
  cfg.sigmaMTDeg,legacy),maxpoolfit.circularTuning(theta,0, ...
  cfg.sigmaMTDeg),'AbsTol',1e-12);
additive=maxpoolfit.configureCoherenceResponse(cfg,"additive");
verifyEqual(testCase,additive.coherenceResponseMode,"additive");
verifyEqual(testCase,additive.modelVersion, ...
  "DenseHardMax-v5-shared-prestep-adaptation");
end

function testPowerAdaptationIsNestedAndDirectional(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureCoherenceResponse(cfg,"additive");
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A,P]=maxpoolfit.deterministicResponse(cache,1,cfg);
Dlinear=maxpoolfit.decisionResponse(B,A,P,0.56,cfg);
verifyEqual(testCase,Dlinear,B-0.56.*A,'AbsTol',1e-12);

cfg=maxpoolfit.configureAdaptationExponent(cfg,1,40);
verifyEqual(testCase,cfg.modelVersion, ...
  "DenseHardMax-v8-power-prestep-adaptation");
Dpower=maxpoolfit.decisionResponse(B,A,P,0.56,cfg);
verifyEqual(testCase,Dpower,B-0.56.*A.*(A./40),'AbsTol',1e-12);
verifyGreaterThan(testCase,max(maxpoolfit.adaptationDrive(A,cfg),[],'all')- ...
  min(maxpoolfit.adaptationDrive(A,cfg),[],'all'),1);
end

function testDirectionalChangeIsBeforeMax(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A,P]=maxpoolfit.deterministicResponse(cache,1:3,cfg);
Dabsolute=maxpoolfit.decisionResponse(B,A,P,0.25,cfg);
verifyEqual(testCase,Dabsolute,B-0.25.*A,'AbsTol',1e-12);
cfg=maxpoolfit.configureDecisionRule(cfg,"directional-change-max");
verifyEqual(testCase,cfg.modelVersion, ...
  "DenseHardMax-v6-directional-change-random-dot-redistribution");
Dchange=maxpoolfit.decisionResponse(B,A,P,0.25,cfg);
verifyEqual(testCase,Dchange,B-P-0.25.*A,'AbsTol',1e-12);
end

function testChoiceProbabilityBounds(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.batchSize=4;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
p=maxpoolfit.choiceProbability(cache,1:5,0.1,8,cfg);
verifySize(testCase,p,[5 1]);
verifyGreaterThan(testCase,min(p),0);
verifyLessThan(testCase,max(p),1);
end

function testAdaptationSuppressesDriftTemplate(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[B,A]=maxpoolfit.deterministicResponse(cache,1,cfg);
D0=B; D1=B-0.2.*A;
[~,i0]=min(abs(cfg.directionGridDeg-0));
verifyLessThan(testCase,D1(i0,1,1),D0(i0,1,1));
end

function testZeroAdaptationRetainsEntireStepResponse(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
t=find(string(cache.trials.domain)=="IDS",1);
[B,A]=maxpoolfit.deterministicResponse(cache,t,cfg);
theta=cfg.directionGridDeg(:);
dirs=double(cache.streamDirDeg{t,1}(:));
cohs=double(cache.stepCohPC{t,1}(:));
preCohs=double(cache.preCohPC{t,1}(:));
expectedR=maxpoolfit.coherenceTuning(theta,dirs,cfg.sigmaMTDeg,cfg)*cohs;
expectedT=maxpoolfit.coherenceTuning(theta,dirs,cfg.sigmaAdaptDeg,cfg)*preCohs;
verifyEqual(testCase,B(:,1,1),expectedR,'AbsTol',1e-12);
verifyEqual(testCase,A(:,1,1),expectedT,'AbsTol',1e-12);
verifyEqual(testCase,B(:,1,1)-0.*A(:,1,1),expectedR,'AbsTol',1e-12);
end

function testProbePedestalsAdaptProbeDirectionMechanisms(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
t=find(string(cache.trials.domain)=="IDQ",1);
cache.preCohPC{t,1}(1)=0;
[B,A]=maxpoolfit.deterministicResponse(cache,t,cfg);
dirs=double(cache.streamDirDeg{t,1}(:));
probeDir=dirs(2);
[~,iProbe]=min(abs(mod(cfg.directionGridDeg-probeDir+180,360)-180));
verifyGreaterThan(testCase,A(iProbe,1,1),0);
verifyLessThan(testCase,B(iProbe,1,1)-0.5*A(iProbe,1,1),B(iProbe,1,1));
end

function testPrestepAdaptationIsDirectionalVector(testCase)
cfg=maxpoolfit.defaultConfig();
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
[~,A]=maxpoolfit.deterministicResponse(cache,1,cfg);
verifyGreaterThan(testCase,max(A(:,1,1))-min(A(:,1,1)),1);
end

function testFitUsesDomainNoiseAndSharedAdaptation(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
fit=maxpoolfit.fitModel(cache,cfg,'NDraws',8);
verifyEqual(testCase,fit.adaptation,cfg.adaptation);
verifyEqual(testCase,fit.domains,["IDR"; "IDQ"; "IDS"]);
verifySize(testCase,fit.noiseSDPC,[3 1]);
verifyEqual(testCase,fit.domainFit.domain,fit.domains);
verifyEqual(testCase,sum(fit.domainFit.nTrials),fit.nTrials);
verifyFalse(testCase,isfield(fit,'animals'));
end

function testFitAcceptsSpecifiedSharedAdaptation(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.adaptation=0.43;
cfg.nNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
fit=maxpoolfit.fitModel(cache,cfg,'NDraws',8);
verifyEqual(testCase,fit.adaptation,0.43,'AbsTol',1e-12);
verifyEqual(testCase,fit.config.adaptation,0.43,'AbsTol',1e-12);
end

function testPredictionSharesNoiseAcrossAnimalsWithinDomain(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
fit=maxpoolfit.fitModel(cache,cfg,'NDraws',8);
prediction=maxpoolfit.predict(cache,fit,cfg,'TrialMask',fit.trialMask,'NDraws',8);
verifyEqual(testCase,prediction.domains,fit.domains);
verifyEqual(testCase,prediction.noiseSDPC,fit.noiseSDPC);
verifyFalse(testCase,isfield(prediction,'animals'));
verifyEqual(testCase,numel(prediction.winnerDetails),3);
end

function testCorrectedAdaptationProfile(testCase)
cfg=maxpoolfit.defaultConfig(); cfg.nNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
P=maxpoolfit.profileAdaptation(cache,cfg,'NDraws',8, ...
  'AdaptationGrid',[0 0.5 1]);
verifyEqual(testCase,P.profile.adaptation,[0;0.5;1]);
verifyEqual(testCase,size(P.profile.noiseSDPC),[3 3]);
verifyEqual(testCase,size(P.profile.domainLogLikelihood),[3 3]);
verifyEqual(testCase,P.domains,["IDR";"IDQ";"IDS"]);
verifyEqual(testCase,P.nTrials,sum(cfg.primarySelector(cache.trials)));
verifyTrue(testCase,all(isfinite(P.profile.logLikelihood)));
end

function testDirectionalChangeProfile(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureDecisionRule(cfg,"directional-change-max");
cfg.nNoiseDraws=8; cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
P=maxpoolfit.profileAdaptation(cache,cfg,'NDraws',8, ...
  'AdaptationGrid',[0 0.5 1]);
verifyEqual(testCase,P.modelVersion, ...
  "DenseHardMax-v6-directional-change-random-dot-redistribution");
verifyEqual(testCase,P.config.decisionRule,"directional-change-max");
verifyTrue(testCase,all(isfinite(P.profile.logLikelihood)));
end


function testAdaptationExponentProfile(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureCoherenceResponse(cfg,"additive"); cfg.nNoiseDraws=8;
cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
P=maxpoolfit.profileAdaptationExponent(cache,cfg,'NDraws',8, ...
  'ExponentGrid',[0 0.5 1]);
verifyEqual(testCase,P.profile.adaptationExponent,[0;0.5;1]);
verifyEqual(testCase,size(P.profile.noiseSDPC),[3 3]);
verifyEqual(testCase,P.baseModelVersion, ...
  "DenseHardMax-v5-shared-prestep-adaptation");
verifyTrue(testCase,all(isfinite(P.profile.logLikelihood)));
verifyEqual(testCase,P.adaptation,cfg.adaptation,'AbsTol',1e-12);
verifyEqual(testCase,P.adaptationReferencePC,40,'AbsTol',1e-12);
end

function testAdaptationGainProfile(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureCoherenceResponse(cfg,"additive"); cfg.nNoiseDraws=8;
cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
P=maxpoolfit.profileAdaptationGain(cache,cfg,'NDraws',8, ...
  'GainGrid',[0 0.5 1]);
verifyEqual(testCase,P.profile.adaptation, ...
  repmat(cfg.adaptation,3,1),'AbsTol',1e-12);
verifyEqual(testCase,P.profile.adaptationGainStrength,[0;0.5;1]);
verifyEqual(testCase,size(P.profile.noiseSDPC),[3 3]);
verifyEqual(testCase,P.baseModelVersion, ...
  "DenseHardMax-v5-shared-prestep-adaptation");
verifyTrue(testCase,all(isfinite(P.profile.logLikelihood)));
verifyEqual(testCase,P.adaptation,cfg.adaptation,'AbsTol',1e-12);
verifyEqual(testCase,P.adaptationReferencePC,40,'AbsTol',1e-12);
end

function testAdaptationGain2DProfile(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureCoherenceResponse(cfg,"additive"); cfg.nNoiseDraws=8;
cfg.batchSize=8; cfg.verbose=false;
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
P=maxpoolfit.profileAdaptationGain(cache,cfg,'NDraws',8, ...
  'Adaptation',0.56,'AdaptationGrid',[0.3 0.56], ...
  'GainGrid',[0 0.5]);
verifyEqual(testCase,height(P.profile),4);
verifyEqual(testCase,unique(P.profile.adaptation),[0.3;0.56], ...
  'AbsTol',1e-12);
verifyEqual(testCase,unique(P.profile.adaptationGainStrength),[0;0.5], ...
  'AbsTol',1e-12);
verifyTrue(testCase,all(isfinite(P.profile.logLikelihood)));
reference=P.profile.adaptation==0.56 & ...
  P.profile.adaptationGainStrength==0;
verifyEqual(testCase,P.profile.deltaLogLFromLinear(reference),0, ...
  'AbsTol',1e-12);
end
