function validatePairedProbeSimulation()
% VALIDATEPAIREDPROBESIMULATION  Fast structural and operating-point checks.

cfg=pairedProbeDefaults();
cfg.OffsetsDeg=[10 120 179];
cfg.BankSpacingDeg=2;
cfg.NTrials=4000;
cfg.NCalibrationTrials=8000;
cfg.NoiseGenerationBatchSize=1000;
cfg.Visible='off';
cfg.SaveOutputs=false;
results=runPairedProbeOffsetSweep(cfg);

assert(height(results.conditionTable)==6,'Expected three offsets by two yoking modes.');
assert(abs(results.noiseCalibration.performance-cfg.TargetPerformance)<0.02, ...
    'The one-time internal-noise calibration missed its performance target.');
assert(all(results.conditionTable.StepCoherence==cfg.FixedStepCoherence), ...
    'Conditions were not analyzed at the common fixed drift increment.');
assert(all(isfinite(results.conditionTable.Threshold75)), ...
    'A fixed-noise 75%% threshold is nonfinite.');
assert(all(results.conditionTable.ProbitConverged),'A probit regression did not converge.');
assert(all(abs(results.conditionTable.ChangedDriftWinner+ ...
    results.conditionTable.ChangedProbeWinner+ ...
    results.conditionTable.ChangedOtherWinner-1)<1e-10), ...
    'Changed-patch winner categories do not sum to one.');
assert(all(isfinite(results.decompositionTable.DecisionTotal)), ...
    'Nonfinite decomposition value detected.');

bank=makeDensePairedTrialBank(200,cfg,cfg.Seed+77);
cfgA=results.config;
modelA=buildPairedProbeResponses(bank,45,false,cfgA);
cfgB=cfgA;
cfgB.BaseDriftCoherence=40;
cfgB.DriftPedestal=13;
cfgB.ProbePedestal=22;
modelB=buildPairedProbeResponses(bank,45,false,cfgB);
assert(isequal(modelA.responseChangedBase,modelB.responseChangedBase) && ...
    isequal(modelA.responseNoChange,modelB.responseNoChange), ...
    'Constant base coherences or pedestals survived mechanism-wise differencing.');
close(results.figureHandles(ishandle(results.figureHandles)));
fprintf('validatePairedProbeSimulation: all checks passed.\n');
end
