function validatePairedProbePoissonSweep()
% VALIDATEPAIREDPROBEPOISSONSWEEP  Fast structural and calibration self-check.

cfg=pairedProbePoissonDefaults();
cfg.BaseDriftCoherences=[7.071 40];
cfg.PairedConfig.OffsetsDeg=[45 120 180];
cfg.PairedConfig.BankSpacingDeg=2;
cfg.PairedConfig.NTrials=2000;
cfg.PairedConfig.NCalibrationTrials=4000;
cfg.PairedConfig.CalibrationTolerance=0.015;
cfg.PairedConfig.CalibrationMaxIterations=22;
cfg.PairedConfig.CreateFigures=false;
cfg.PairedConfig.SaveOutputs=false;
cfg.IncludeHomoscedasticReference=false;
cfg.CreateFigures=false;
cfg.SaveOutputs=false;
cfg.Visible='off';
results=runPairedProbePoissonSweep(cfg);
T=results.conditionTable;
assert(height(T)==12,'Unexpected number of conditions.');
assert(results.poolCalibration.effectivePoolSize>0 && ...
    isfinite(results.poolCalibration.effectivePoolSize),'Invalid pool size.');
assert(abs(results.poolCalibration.performance- ...
    cfg.PairedConfig.TargetPerformance)<0.04,'Pool calibration missed target.');
assert(all(isfinite(T.Threshold75)) && all(T.Threshold75>0), ...
    'Invalid 75%% threshold.');
assert(all(T.EffectivePoolSize==T.EffectivePoolSize(1)), ...
    'Effective pool size was not held fixed.');
assert(all(T.ProbitConverged),'At least one probit fit did not converge.');
assert(max([T.PreClampFraction;T.ChangedClampFraction;T.NoChangeClampFraction])<0.02, ...
    'Rate positivity bound affects too many samples.');
fprintf('Poisson paired-probe validation passed. Calibrated N_eff = %.4g.\n', ...
    results.poolCalibration.effectivePoolSize);
end
