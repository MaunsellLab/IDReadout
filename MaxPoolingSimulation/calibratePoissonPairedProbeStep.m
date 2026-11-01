function step = calibratePoissonPairedProbeStep(bank,offsetDeg,yoked, ...
    baseDriftCoherence,effectivePoolSize,cfg)
% CALIBRATEPOISSONPAIREDPROBESTEP  75% threshold at a fixed effective pool.

pcfg=cfg.PairedConfig;
performanceAt=@performanceForStep;
lower=0;
upper=pcfg.CalibrationInitialUpper;
while performanceAt(upper)<pcfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket target performance.'); end
end
for iter=1:pcfg.CalibrationMaxIterations
    mid=(lower+upper)/2;
    p=performanceAt(mid);
    if abs(p-pcfg.TargetPerformance)<=pcfg.CalibrationTolerance
        step=mid;
        return;
    elseif p<pcfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end
step=(lower+upper)/2;

    function p=performanceForStep(candidateStep)
        model=preparePoissonPairedProbeCondition(bank,offsetDeg,yoked, ...
            candidateStep,baseDriftCoherence,cfg);
        scale=single(1/sqrt(effectivePoolSize));
        changed=max(model.meanChanged+scale*model.noiseChangedUnit,[],2);
        noChange=max(model.meanNoChange+scale*model.noiseNoChangeUnit,[],2);
        p=mean(changed>noChange);
    end
end
