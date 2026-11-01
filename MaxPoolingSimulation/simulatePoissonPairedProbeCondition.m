function sim = simulatePoissonPairedProbeCondition(bank,offsetDeg,yoked, ...
    stepCoherence,baseDriftCoherence,effectivePoolSize,cfg)
% SIMULATEPOISSONPAIREDPROBECONDITION  Choices and gains with count variance.

model=preparePoissonPairedProbeCondition(bank,offsetDeg,yoked, ...
    stepCoherence,baseDriftCoherence,cfg);
scale=single(1/sqrt(effectivePoolSize));
responseChanged=model.meanChanged+scale*model.noiseChangedUnit;
responseNoChange=model.meanNoChange+scale*model.noiseNoChangeUnit;
[valueChanged,winnerChanged]=max(responseChanged,[],2);
[valueNoChange,winnerNoChange]=max(responseNoChange,[],2);
correct=valueChanged>valueNoChange;
fit=fitProbitIRLS(model.predictors,correct);

if yoked
    cTotal=[0;1]; cPer=[0;0.5]; yokingName='Yoked';
else
    cTotal=[0;1;1]; cPer=[0;0.5;0.5]; yokingName='Independent';
end
covariance=fit.covariance;
gain=struct();
gain.drift=fit.beta(1);
gain.driftSE=sqrt(max(covariance(1,1),0));
gain.probeTotal=cTotal'*fit.beta;
gain.probeTotalSE=sqrt(max(cTotal'*covariance*cTotal,0));
gain.probePerStream=cPer'*fit.beta;
gain.probePerStreamSE=sqrt(max(cPer'*covariance*cPer,0));

winners=classifyDenseWinners(winnerChanged,winnerNoChange, ...
    bank.directionsDeg,model.streamDirections, ...
    cfg.PairedConfig.WinnerAssociationRadiusDeg);

sim=struct();
sim.offsetDeg=offsetDeg;
sim.yoked=yoked;
sim.yoking=yokingName;
sim.baseDriftCoherence=baseDriftCoherence;
sim.stepCoherence=stepCoherence;
sim.effectivePoolSize=effectivePoolSize;
sim.performance=mean(correct);
sim.fit=fit;
sim.gain=gain;
sim.winners=winners;
sim.rateDiagnostics=rmfield(model,{'W','meanChanged','meanNoChange', ...
    'noiseChangedUnit','noiseNoChangeUnit','predictors', ...
    'streamDirections','offsetDeg','yoked'});
end
