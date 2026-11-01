function threshold = calibrateSplitSignalModelThreshold(bank,modelName,allocation,cfg)
% CALIBRATESPLITSIGNALMODELTHRESHOLD  Bisection to target performance.

bounds=cfg.ThresholdInitialBounds;
low=bounds(1); high=bounds(2);
pLow=performanceAt(low); pHigh=performanceAt(high);
while pHigh<cfg.TargetPerformance && high<cfg.MaximumTotalStepCoherence
    high=min(cfg.MaximumTotalStepCoherence,1.2*high); pHigh=performanceAt(high);
end
if pLow>=cfg.TargetPerformance
    threshold=low; return;
end
if pHigh<cfg.TargetPerformance
    threshold=NaN; return;
end
for iter=1:cfg.ThresholdMaxIterations
    middle=0.5*(low+high); pMiddle=performanceAt(middle);
    if pMiddle>=cfg.TargetPerformance, high=middle; else, low=middle; end
    if high-low<=cfg.ThresholdTolerance, break; end
end
threshold=0.5*(low+high);

    function p=performanceAt(step)
        model=buildSplitSignalDesign(bank,cfg.SeparationDeg,allocation,step, ...
            cfg.StepNoiseHalfAmplitude,cfg.PrestepPedestalFraction, ...
            cfg.BalancePrestep,cfg.AdaptationFraction,cfg);
        choice=simulateSplitSignalPoolingChoice(model,bank,modelName, ...
            cfg.TrueInternalNoiseSD);
        p=mean(choice);
    end
end
