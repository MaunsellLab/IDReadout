function sim = simulateTemporalMaxCondition(bank,offsetDeg,yoked, ...
    stepCoherence,baseDriftCoherence,internalStepNoiseSD,cfg)
% SIMULATETEMPORALMAXCONDITION  Temporal evidence, gains, and winner dynamics.

model=buildTemporalMaxResponses(bank,offsetDeg,yoked,stepCoherence, ...
    baseDriftCoherence,cfg);
scale=single(internalStepNoiseSD);
preChanged=bsxfun(@plus,model.preMean,scale*model.preChangedUnit);
preNoChange=bsxfun(@plus,model.preMean,scale*model.preNoChangeUnit);
stepChanged=model.stepChangedMean+scale*model.stepChangedUnit;
stepNoChange=model.stepNoChangeMean+scale*model.stepNoChangeUnit;
[preValueC,preWinnerC]=max(preChanged,[],2);
[preValueN,preWinnerN]=max(preNoChange,[],2);
[stepValueC,stepWinnerC]=max(stepChanged,[],2);
[stepValueN,stepWinnerN]=max(stepNoChange,[],2);
evidenceChanged=stepValueC-preValueC;
evidenceNoChange=stepValueN-preValueN;
correct=evidenceChanged>evidenceNoChange;
fit=fitProbitIRLS(model.predictors,correct);

if yoked
    cTotal=[0;1]; cPer=[0;0.5]; yokingName='Yoked';
else
    cTotal=[0;1;1]; cPer=[0;0.5;0.5]; yokingName='Independent';
end
covariance=fit.covariance;
gain=struct('drift',fit.beta(1), ...
    'driftSE',sqrt(max(covariance(1,1),0)), ...
    'probeTotal',cTotal'*fit.beta, ...
    'probeTotalSE',sqrt(max(cTotal'*covariance*cTotal,0)), ...
    'probePerStream',cPer'*fit.beta, ...
    'probePerStreamSE',sqrt(max(cPer'*covariance*cPer,0)));

winners=classifyDenseWinners(stepWinnerC,stepWinnerN, ...
    bank.directionsDeg,model.streamDirections, ...
    cfg.PairedConfig.WinnerAssociationRadiusDeg);
transitions=struct();
transitions.changed=transitionMatrix(preWinnerC,stepWinnerC,bank.directionsDeg, ...
    model.streamDirections,cfg.PairedConfig.WinnerAssociationRadiusDeg);
transitions.noChange=transitionMatrix(preWinnerN,stepWinnerN,bank.directionsDeg, ...
    model.streamDirections,cfg.PairedConfig.WinnerAssociationRadiusDeg);
transitions.changedSame=mean(preWinnerC==stepWinnerC);
transitions.noChangeSame=mean(preWinnerN==stepWinnerN);

winnerDirection=computeWinnerDirectionDiagnostic(stepWinnerC,stepWinnerN, ...
    bank.directionsDeg,model.W,cfg.WinnerDirectionBinWidthDeg);

driftDecomp=decomposeTemporalMaxInfluence(preWinnerC,preWinnerN, ...
    stepChanged,stepNoChange,model.W(1,:),cfg.PairedConfig.CounterfactualPerturbation);
if yoked
    probeDecomp=decomposeTemporalMaxInfluence(preWinnerC,preWinnerN, ...
        stepChanged,stepNoChange,model.W(2,:)+model.W(3,:), ...
        cfg.PairedConfig.CounterfactualPerturbation);
    probeDecomp=scaleDecomp(probeDecomp,0.5);
else
    p1=decomposeTemporalMaxInfluence(preWinnerC,preWinnerN, ...
        stepChanged,stepNoChange,model.W(2,:), ...
        cfg.PairedConfig.CounterfactualPerturbation);
    p2=decomposeTemporalMaxInfluence(preWinnerC,preWinnerN, ...
        stepChanged,stepNoChange,model.W(3,:), ...
        cfg.PairedConfig.CounterfactualPerturbation);
    probeDecomp=averageDecomp(p1,p2);
end

sim=struct('offsetDeg',offsetDeg,'yoked',yoked,'yoking',yokingName, ...
    'baseDriftCoherence',baseDriftCoherence,'stepCoherence',stepCoherence, ...
    'performance',mean(correct),'fit',fit,'gain',gain,'winners',winners, ...
    'transitions',transitions,'winnerDirection',winnerDirection);
sim.decomposition.drift=driftDecomp;
sim.decomposition.probePerStream=probeDecomp;
end

function matrix=transitionMatrix(preIndex,stepIndex,directions,streams,radius)
preCategory=winnerCategory(preIndex,directions,streams,radius);
stepCategory=winnerCategory(stepIndex,directions,streams,radius);
matrix=accumarray([preCategory,stepCategory],1,[3 3])/numel(preCategory);
end

function category=winnerCategory(index,directions,streams,radius)
winnerDirection=directions(index(:));
driftDistance=abs(mod(winnerDirection-streams(1)+180,360)-180);
probeDistance=min(abs(mod(winnerDirection-streams(2)+180,360)-180), ...
    abs(mod(winnerDirection-streams(3)+180,360)-180));
isDrift=driftDistance<=radius & driftDistance<=probeDistance;
isProbe=~isDrift & probeDistance<=radius;
category=3*ones(numel(index),1); category(isProbe)=2; category(isDrift)=1;
end

function d=scaleDecomp(d,scale)
names=fieldnames(d);
for k=1:numel(names)
    if ~contains(names{k},'Probability'), d.(names{k})=scale*d.(names{k}); end
end
end

function d=averageDecomp(a,b)
d=a; names=fieldnames(a);
for k=1:numel(names), d.(names{k})=0.5*(a.(names{k})+b.(names{k})); end
end
