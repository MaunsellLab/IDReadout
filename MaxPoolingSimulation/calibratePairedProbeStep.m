function step = calibratePairedProbeStep(bank,offsetDeg,yoked,cfg)
% CALIBRATEPAIREDPROBESTEP  Measure the 75% threshold with noise held fixed.
%
% This function does not set an operating point for the analysis.  The single
% internal-noise calibration has already occurred, and the same fixed noise is
% used here at every offset and yoking condition.

model=buildPairedProbeResponses(bank,offsetDeg,yoked,cfg);
driftProfile=model.W(1,:);
noChangeValue=max(model.responseNoChange,[],2);
performanceAt=@(s)mean(max(bsxfun(@plus,model.responseChangedBase, ...
    single(s)*driftProfile),[],2)>noChangeValue);
lower=0; upper=cfg.CalibrationInitialUpper;
while performanceAt(upper)<cfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket target performance.'); end
end
for iter=1:cfg.CalibrationMaxIterations
    mid=(lower+upper)/2; p=performanceAt(mid);
    if abs(p-cfg.TargetPerformance)<=cfg.CalibrationTolerance
        step=mid; return;
    elseif p<cfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end
step=(lower+upper)/2;
end
