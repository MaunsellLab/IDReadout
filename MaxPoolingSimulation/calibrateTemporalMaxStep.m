function step = calibrateTemporalMaxStep(bank,offsetDeg,yoked, ...
    baseDriftCoherence,internalStepNoiseSD,cfg)
% CALIBRATETEMPORALMAXSTEP  Measure the 75% threshold with noise held fixed.

pcfg=cfg.PairedConfig;
model=buildTemporalMaxResponses(bank,offsetDeg,yoked,0,baseDriftCoherence,cfg);
scale=single(internalStepNoiseSD);
preC=max(bsxfun(@plus,model.preMean,scale*model.preChangedUnit),[],2);
preN=max(bsxfun(@plus,model.preMean,scale*model.preNoChangeUnit),[],2);
stepN=max(model.stepNoChangeMean+scale*model.stepNoChangeUnit,[],2);
performanceAt=@(s)mean((max(bsxfun(@plus, ...
    model.stepChangedMean+scale*model.stepChangedUnit,single(s)*model.W(1,:)), ...
    [],2)-preC)>(stepN-preN));
lower=0; upper=pcfg.CalibrationInitialUpper;
while performanceAt(upper)<pcfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket target performance.'); end
end
for iter=1:pcfg.CalibrationMaxIterations
    mid=(lower+upper)/2; p=performanceAt(mid);
    if abs(p-pcfg.TargetPerformance)<=pcfg.CalibrationTolerance
        step=mid; return;
    elseif p<pcfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end
step=(lower+upper)/2;
end
