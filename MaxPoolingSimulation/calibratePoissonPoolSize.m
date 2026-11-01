function calibration = calibratePoissonPoolSize(bank,cfg)
% CALIBRATEPOISSONPOOLSIZE  Fit one effective pool size at the reference base.

pcfg=cfg.PairedConfig;
model=preparePoissonPairedProbeCondition(bank, ...
    pcfg.NoiseCalibrationOffsetDeg,logical(pcfg.NoiseCalibrationYoked), ...
    pcfg.FixedStepCoherence,cfg.ReferenceBaseDriftCoherence,cfg);
performanceAt=@performanceForPool;

lower=1e-4;
pLower=performanceAt(lower);
if pLower>pcfg.TargetPerformance+pcfg.CalibrationTolerance
    error('Even a very small effective pool exceeds target performance.');
end
upper=max(cfg.PoolCalibrationInitialUpper,1);
pUpper=performanceAt(upper);
while pUpper<pcfg.TargetPerformance
    upper=2*upper;
    if upper>cfg.PoolCalibrationMaximum
        error('Could not bracket the effective pool size.');
    end
    pUpper=performanceAt(upper);
end

mid=upper; pMid=pUpper;
for iter=1:pcfg.CalibrationMaxIterations
    mid=sqrt(lower*upper);
    pMid=performanceAt(mid);
    if abs(pMid-pcfg.TargetPerformance)<=pcfg.CalibrationTolerance
        break;
    elseif pMid<pcfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end

calibration=struct();
calibration.effectivePoolSize=mid;
calibration.performance=pMid;
calibration.performanceAtTinyPool=pLower;
calibration.performanceWithoutSpikeNoise=performanceAt(Inf);
calibration.baseDriftCoherence=cfg.ReferenceBaseDriftCoherence;
calibration.offsetDeg=pcfg.NoiseCalibrationOffsetDeg;
calibration.yoked=logical(pcfg.NoiseCalibrationYoked);
calibration.stepCoherence=pcfg.FixedStepCoherence;
calibration.randomMotionRateSpPerSec=cfg.RandomMotionRateSpPerSec;
calibration.alignedRateSlopeSpPerSecPerPct= ...
    cfg.AlignedRateSlopeSpPerSecPerPct;
calibration.preStepDurationSec=cfg.PreStepDurationSec;
calibration.stepDurationSec=cfg.StepDurationSec;
calibration.fanoFactor=cfg.FanoFactor;
calibration.correlationLengthDeg=pcfg.InternalCorrelationLengthDeg;
calibration.effectiveNoiseRank=circularNoiseEffectiveRank( ...
    numel(bank.directionsDeg),pcfg);
calibration.preRateDrift=model.preRateDrift;
calibration.preRateProbe=model.preRateProbe;

    function p=performanceForPool(poolSize)
        if isinf(poolSize)
            changed=max(model.meanChanged,[],2);
            noChange=max(model.meanNoChange,[],2);
        else
            scale=single(1/sqrt(poolSize));
            changed=max(model.meanChanged+scale*model.noiseChangedUnit,[],2);
            noChange=max(model.meanNoChange+scale*model.noiseNoChangeUnit,[],2);
        end
        p=mean(changed>noChange);
    end
end
