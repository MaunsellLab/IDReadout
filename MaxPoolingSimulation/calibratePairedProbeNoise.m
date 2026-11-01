function calibration = calibratePairedProbeNoise(bank,cfg)
% CALIBRATEPAIREDPROBENOISE  Set one internal-noise SD at a fixed 20% step.
%
% External stream-change noise remains fixed in coherence units.  This
% routine varies only the additive, direction-correlated mechanism-noise SD
% until FixedStepCoherence produces TargetPerformance at the single reference
% offset/yoking condition.  The returned SD is then held fixed for every
% offset and both yoking conditions.

offsetDeg=cfg.NoiseCalibrationOffsetDeg;
yoked=logical(cfg.NoiseCalibrationYoked);
step=cfg.FixedStepCoherence;

zeroNoiseCfg=cfg;
zeroNoiseCfg.InternalNoiseSD=0;
model=buildPairedProbeResponses(bank,offsetDeg,yoked,zeroNoiseCfg);
changedSensory=bsxfun(@plus,model.responseChangedBase, ...
    single(step)*model.W(1,:));
noChangeSensory=model.responseNoChange;

performanceAt=@performanceForSD;
pZero=performanceAt(0);
if pZero<cfg.TargetPerformance-cfg.CalibrationTolerance
    error(['External stream noise alone puts performance below the target. ' ...
        'An additive internal-noise calibration cannot correct this.']);
end

lower=0;
upper=max(cfg.InternalNoiseSD,cfg.CalibrationInitialUpper);
pUpper=performanceAt(upper);
while pUpper>cfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket the target internal-noise SD.'); end
    pUpper=performanceAt(upper);
end

mid=upper;
pMid=pUpper;
for iter=1:cfg.CalibrationMaxIterations
    mid=(lower+upper)/2;
    pMid=performanceAt(mid);
    if abs(pMid-cfg.TargetPerformance)<=cfg.CalibrationTolerance
        break;
    elseif pMid>cfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end

calibration=struct();
calibration.internalNoiseSD=mid;
calibration.performance=pMid;
calibration.performanceWithNoInternalNoise=pZero;
calibration.offsetDeg=offsetDeg;
calibration.yoked=yoked;
calibration.stepCoherence=step;

    function p=performanceForSD(noiseSD)
        changed=max(changedSensory+single(noiseSD)*bank.internalChangedUnit,[],2);
        noChange=max(noChangeSensory+single(noiseSD)*bank.internalNoChangeUnit,[],2);
        p=mean(changed>noChange);
    end
end
