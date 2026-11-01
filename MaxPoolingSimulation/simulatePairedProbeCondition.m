function sim = simulatePairedProbeCondition(bank,offsetDeg,yoked,step,cfg)
% SIMULATEPAIREDPROBECONDITION  Choices, gains, winners, and decompositions.

model=buildPairedProbeResponses(bank,offsetDeg,yoked,cfg);
responseChanged=bsxfun(@plus,model.responseChangedBase,single(step)*model.W(1,:));
[valueChanged,winnerChanged]=max(responseChanged,[],2);
[valueNoChange,winnerNoChange]=max(model.responseNoChange,[],2);
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
    bank.directionsDeg,model.streamDirections,cfg.WinnerAssociationRadiusDeg);

driftDecomp=decomposeStreamInfluence(responseChanged,model.responseNoChange, ...
    model.W(1,:),cfg.CounterfactualPerturbation);
if yoked
    probePair=decomposeStreamInfluence(responseChanged,model.responseNoChange, ...
        model.W(2,:)+model.W(3,:),cfg.CounterfactualPerturbation);
    probePer=scaleDecomposition(probePair,0.5);
else
    probe1=decomposeStreamInfluence(responseChanged,model.responseNoChange, ...
        model.W(2,:),cfg.CounterfactualPerturbation);
    probe2=decomposeStreamInfluence(responseChanged,model.responseNoChange, ...
        model.W(3,:),cfg.CounterfactualPerturbation);
    probePer=averageDecompositions(probe1,probe2);
end

sim=struct();
sim.offsetDeg=offsetDeg;
sim.yoked=yoked;
sim.yoking=yokingName;
sim.stepCoherence=step;
sim.performance=mean(correct);
sim.fit=fit;
sim.gain=gain;
sim.winners=winners;
sim.decomposition.drift=driftDecomp;
sim.decomposition.probePerStream=probePer;
end

function d=scaleDecomposition(d,scale)
fields={'changedFixed','changedSwitch','noChangeFixed','noChangeSwitch', ...
    'decisionFixed','decisionSwitch','decisionTotal'};
for k=1:numel(fields), d.(fields{k})=scale*d.(fields{k}); end
end

function d=averageDecompositions(a,b)
d=a;
fields={'changedFixed','changedSwitch','noChangeFixed','noChangeSwitch', ...
    'decisionFixed','decisionSwitch','decisionTotal', ...
    'changedSwitchProbability','noChangeSwitchProbability'};
for k=1:numel(fields), d.(fields{k})=0.5*(a.(fields{k})+b.(fields{k})); end
end

