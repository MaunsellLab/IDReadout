function calibration = calibrateTemporalMaxNoise(bank,cfg)
% CALIBRATETEMPORALMAXNOISE  Fit one step-epoch SD at the reference condition.

pcfg=cfg.PairedConfig;
model=buildTemporalMaxResponses(bank,pcfg.NoiseCalibrationOffsetDeg, ...
    logical(pcfg.NoiseCalibrationYoked),pcfg.FixedStepCoherence, ...
    cfg.ReferenceBaseDriftCoherence,cfg);
performanceAt=@performanceForSD;
pZero=performanceAt(0);
if pZero<pcfg.TargetPerformance-pcfg.CalibrationTolerance
    error('External stream noise alone puts performance below target.');
end
lower=0; upper=max(cfg.InternalStepNoiseSD,pcfg.CalibrationInitialUpper);
while performanceAt(upper)>pcfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket target internal-noise SD.'); end
end
mid=upper; pMid=performanceAt(mid);
for iter=1:pcfg.CalibrationMaxIterations
    mid=(lower+upper)/2; pMid=performanceAt(mid);
    if abs(pMid-pcfg.TargetPerformance)<=pcfg.CalibrationTolerance
        break;
    elseif pMid>pcfg.TargetPerformance
        lower=mid;
    else
        upper=mid;
    end
end
calibration=struct('internalStepNoiseSD',mid,'performance',pMid, ...
    'performanceWithNoInternalNoise',pZero, ...
    'internalPreStepNoiseSD',mid*sqrt(cfg.StepDurationSec/cfg.PreStepDurationSec), ...
    'baseDriftCoherence',cfg.ReferenceBaseDriftCoherence, ...
    'offsetDeg',pcfg.NoiseCalibrationOffsetDeg, ...
    'yoked',logical(pcfg.NoiseCalibrationYoked), ...
    'stepCoherence',pcfg.FixedStepCoherence, ...
    'preStepDurationSec',cfg.PreStepDurationSec, ...
    'stepDurationSec',cfg.StepDurationSec, ...
    'correlationLengthDeg',pcfg.InternalCorrelationLengthDeg, ...
    'effectiveNoiseRank',circularNoiseEffectiveRank(numel(bank.directionsDeg),pcfg));

    function p=performanceForSD(sd)
        scale=single(sd);
        preC=max(bsxfun(@plus,model.preMean,scale*model.preChangedUnit),[],2);
        stepC=max(model.stepChangedMean+scale*model.stepChangedUnit,[],2);
        preN=max(bsxfun(@plus,model.preMean,scale*model.preNoChangeUnit),[],2);
        stepN=max(model.stepNoChangeMean+scale*model.stepNoChangeUnit,[],2);
        p=mean((stepC-preC)>(stepN-preN));
    end
end
