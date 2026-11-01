function calibration = calibrateAdaptationNoise(bank,cfg)
% CALIBRATEADAPTATIONNOISE  One fixed internal-noise SD at the reference state.

pcfg=cfg.PairedConfig;
offsetDeg=pcfg.NoiseCalibrationOffsetDeg;
yoked=logical(pcfg.NoiseCalibrationYoked);
base=cfg.ReferenceBaseDriftCoherence;
aftereffect=cfg.ReferenceAftereffectAt60;
step=pcfg.FixedStepCoherence;

model=buildAdaptationResponses(bank,offsetDeg,yoked,base,aftereffect,0,cfg);
changedSensory=bsxfun(@plus,model.responseChangedBase, ...
    single(step)*model.W(1,:));
noChangeSensory=model.responseNoChange;
performanceAt=@performanceForSD;
pZero=performanceAt(0);
if pZero<pcfg.TargetPerformance-pcfg.CalibrationTolerance
    error(['External stream noise alone puts performance below the target. ' ...
        'An additive internal-noise calibration cannot correct this.']);
end

lower=0; upper=max(pcfg.InternalNoiseSD,pcfg.CalibrationInitialUpper);
pUpper=performanceAt(upper);
while pUpper>pcfg.TargetPerformance
    upper=2*upper;
    if upper>1e6, error('Could not bracket the target internal-noise SD.'); end
    pUpper=performanceAt(upper);
end
mid=upper; pMid=pUpper;
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
calibration=struct('internalNoiseSD',mid,'performance',pMid, ...
    'performanceWithNoInternalNoise',pZero,'baseDriftCoherence',base, ...
    'aftereffectAt60',aftereffect, ...
    'actualAftereffect',model.actualAftereffect,'offsetDeg',offsetDeg, ...
    'yoked',yoked,'stepCoherence',step);

    function p=performanceForSD(noiseSD)
        changed=max(changedSensory+single(noiseSD)*bank.internalChangedUnit,[],2);
        noChange=max(noChangeSensory+single(noiseSD)*bank.internalNoChangeUnit,[],2);
        p=mean(changed>noChange);
    end
end
