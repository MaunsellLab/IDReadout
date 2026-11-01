function step = calibrateAdaptationStep(bank,offsetDeg,yoked, ...
    baseDriftCoherence,aftereffectAt60,internalNoiseSD,cfg)
% CALIBRATEADAPTATIONSTEP  Measure a 75% threshold with fixed internal noise.

pcfg=cfg.PairedConfig;
model=buildAdaptationResponses(bank,offsetDeg,yoked,baseDriftCoherence, ...
    aftereffectAt60,internalNoiseSD,cfg);
driftProfile=model.W(1,:);
noChangeValue=max(model.responseNoChange,[],2);
performanceAt=@(s)mean(max(bsxfun(@plus,model.responseChangedBase, ...
    single(s)*driftProfile),[],2)>noChangeValue);
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
