function tests=testIndependentNoiseProfile
tests=functiontests(localfunctions);
end

function testZeroFractionExactlyNestsPriorCore(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
cfg=maxpoolfit.defaultConfig(); cfg.batchSize=8;
idx=(1:8).'; [B,A,P]=maxpoolfit.deterministicResponse(cache,idx,cfg);
z1=maxpoolfit.noiseDraws(cfg,8,0); z2=maxpoolfit.noiseDraws(cfg,8,1000003);
[oldP,oldD]=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,P,.4,10,cfg,z1,z2,true);
newCfg=maxpoolfit.configureNoiseMixture(cfg,0,2);
[newP,newD]=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,P,.4,10,newCfg,z1,z2,true);
verifyEqual(testCase,newP,oldP,'AbsTol',0);
verifyEqual(testCase,newD.winnerMechanismIndexRDK1, ...
  oldD.winnerMechanismIndexRDK1);
verifyEqual(testCase,newD.winnerMechanismIndexRDK2, ...
  oldD.winnerMechanismIndexRDK2);
end

function testResponseVarianceExponentContract(testCase)
cfg=maxpoolfit.defaultConfig(); B=zeros(3,2,1);
B(:,1,1)=[-5;0;10];
cfg=maxpoolfit.configureNoiseMixture(cfg,.2,0);
[s0,F,F0]=maxpoolfit.responseNoiseScale(B,cfg);
verifyEqual(testCase,s0,ones(size(B)));
cfg=maxpoolfit.configureNoiseMixture(cfg,.2,1);
s1=maxpoolfit.responseNoiseScale(B,cfg);
verifyEqual(testCase,s1(:,1,1).^2,F(:,1,1)./F0,'AbsTol',1e-12);
verifyGreaterThan(testCase,s1(3,1,1),s1(2,1,1));
verifyGreaterThan(testCase,s1(2,1,1),s1(1,1,1));
end

function testNegativeLinearProxyIsRectified(testCase)
cfg=maxpoolfit.defaultConfig();
cfg=maxpoolfit.configureNoiseMixture(cfg,.25,1);
F0=100*maxpoolfit.randomDirectionMean(cfg.sigmaMTDeg);
B=zeros(2,2,1); B(1,1,1)=-2*F0;
[scale,fullResponse]=maxpoolfit.responseNoiseScale(B,cfg);
verifyGreaterThan(testCase,fullResponse(1,1,1),0);
verifyTrue(testCase,all(isfinite(scale),'all'));
verifyGreaterThan(testCase,scale(1,1,1),0);
verifyLessThan(testCase,scale(1,1,1),scale(2,1,1));
end

function testLikelihoodOnlyProfileContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
trialGroup=string(cache.trials.domain);
groups=unique(trialGroup,'stable');
adaptationByGroup=table(groups,repmat(.4,numel(groups),1), ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
baseFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
P=maxpoolfit.profileIndependentNoiseMixture(cache,baseFit, ...
  'IndependentFractionGrid',[0 .25], ...
  'VarianceExponentGrid',[0 1],'NDraws',8, ...
  'ComputeGainDiagnostics',false,'ComputeWinnerDiagnostics',false);
verifyEqual(testCase,height(P.overallProfile),4);
verifyEqual(testCase,height(P.groupProfile),12);
verifyEqual(testCase,P.independentFractionGrid,[0;.25]);
verifyEqual(testCase,P.varianceExponentGrid,[0;1]);
verifyTrue(testCase,all(isfinite(P.overallProfile.logLikelihood)));
q=P.overallProfile.independentNoiseFraction==0;
verifyEqual(testCase,P.overallProfile.logLikelihood(q), ...
  repmat(baseFit.logLikelihood,2,1),'AbsTol',1e-10);
end

function testWinnerSummaryContract(testCase)
cache=maxpoolfit.syntheticCache('TrialsPerCondition',1);
trialGroup=string(cache.trials.domain);
groups=unique(trialGroup,'stable');
adaptationByGroup=table(groups,repmat(.4,numel(groups),1), ...
  'VariableNames',{'analysisGroup','adaptation'});
cfg=maxpoolfit.defaultConfig(); cfg.verbose=false; cfg.batchSize=8;
fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'NDraws',8);
prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg, ...
  'NDraws',8,'ReturnWinners',true);
W=maxpoolfit.winnerFractionDiagnostics(cache,prediction,cfg,trialGroup);
verifyFalse(testCase,isempty(W));
verifyTrue(testCase,all(W.changedProbeWinnerFraction>=0 & ...
  W.changedProbeWinnerFraction<=1));
verifyEqual(testCase,W.changedProbeWinnerFraction+ ...
  W.changedDriftWinnerFraction,ones(height(W),1),'AbsTol',1e-12);
verifyEqual(testCase,unique(W.winnerAssignment), ...
  "nearest presented stream direction");
end
