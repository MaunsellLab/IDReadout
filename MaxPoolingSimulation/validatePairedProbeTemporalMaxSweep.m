function validatePairedProbeTemporalMaxSweep()
% VALIDATEPAIREDPROBETEMPORALMAXSWEEP  Fast calibration and table self-check.

cfg=temporalMaxDefaults();
cfg.PairedConfig.OffsetsDeg=[45 120 180];
cfg.PairedConfig.BankSpacingDeg=2;
cfg.PairedConfig.NTrials=2000;
cfg.PairedConfig.NCalibrationTrials=4000;
cfg.PairedConfig.CalibrationTolerance=0.015;
cfg.PairedConfig.CalibrationMaxIterations=22;
cfg.CreateFigures=false; cfg.SaveOutputs=false; cfg.Visible='off';
results=runPairedProbeTemporalMaxSweep(cfg);
T=results.conditionTable; X=results.transitionTable;
W=results.winnerDirectionTable;
assert(height(T)==12,'Unexpected number of conditions.');
assert(height(X)==216,'Unexpected number of transition rows.');
assert(height(W)==12*(360/cfg.WinnerDirectionBinWidthDeg), ...
    'Unexpected number of winner-direction rows.');
assert(abs(results.noiseCalibration.performance- ...
    cfg.PairedConfig.TargetPerformance)<0.04,'Noise calibration missed target.');
assert(all(isfinite(T.Threshold75)) && all(T.Threshold75>0), ...
    'Invalid threshold.');
assert(all(T.ProbitConverged),'At least one probit fit did not converge.');
groups=findgroups(X.BaseDriftCoherence,X.OffsetDeg,X.Yoking,X.Patch);
sums=splitapply(@sum,X.Fraction,groups);
assert(max(abs(sums-1))<1e-10,'Winner transition matrices do not sum to one.');
winnerGroups=findgroups(W.BaseDriftCoherence,W.OffsetDeg,W.Yoking);
winnerSums=splitapply(@sum,W.PatchAverageWinnerFraction,winnerGroups);
assert(max(abs(winnerSums-1))<1e-10, ...
    'Winner-direction distributions do not sum to one.');
fprintf('Temporal-max validation passed. Step-noise SD = %.4g.\n', ...
    results.noiseCalibration.internalStepNoiseSD);
end
