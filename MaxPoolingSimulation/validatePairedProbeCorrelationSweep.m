function validatePairedProbeCorrelationSweep()
% VALIDATEPAIREDPROBECORRELATIONSWEEP  Fast structural sweep checks.

cfg=pairedProbeCorrelationSweepDefaults();
cfg.CorrelationLengthsDeg=[20 60];
cfg.PairedConfig.OffsetsDeg=[45 120 180];
cfg.PairedConfig.BankSpacingDeg=2;
cfg.PairedConfig.NTrials=2500;
cfg.PairedConfig.NCalibrationTrials=5000;
cfg.PairedConfig.NoiseGenerationBatchSize=1000;
cfg.Visible='off';
cfg.SaveOutputs=false;
results=runPairedProbeCorrelationSweep(cfg);

assert(height(results.conditionTable)==12, ...
    'Expected two lengths by three offsets by two yoking modes.');
assert(height(results.noiseCalibrationTable)==2, ...
    'Expected one calibration row per correlation length.');
assert(all(abs(results.noiseCalibrationTable.performance- ...
    cfg.PairedConfig.TargetPerformance)<0.03), ...
    'A one-time noise calibration missed its target.');
rows=sortrows(results.noiseCalibrationTable,'correlationLengthDeg');
assert(rows.effectiveNoiseRank(1)>rows.effectiveNoiseRank(2), ...
    'Longer correlation length did not reduce effective dimensionality.');
assert(all(results.conditionTable.StepCoherence== ...
    cfg.PairedConfig.FixedStepCoherence), ...
    'A sweep condition did not use the common fixed increment.');
assert(all(results.conditionTable.ProbitConverged), ...
    'A probit regression did not converge.');
close(results.figureHandles(ishandle(results.figureHandles)));
fprintf('validatePairedProbeCorrelationSweep: all checks passed.\n');
end
