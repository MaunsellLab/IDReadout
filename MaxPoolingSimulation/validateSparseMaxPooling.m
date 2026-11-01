function validateSparseMaxPooling()
% VALIDATESPARSEMAXPOOLING  Fast checks of alignment, centering, and calibration.
%
% This is a small smoke test, not a statistical validation of the model.

cfg = maxPoolingDefaults();
cfg.NValues = [2 3 6];
cfg.ReferenceN = 3;
cfg.DiagnosticN = 3;
cfg.NTrials = 6000;
cfg.NCalibrationTrials = 12000;
cfg.Visible = 'off';
cfg.SaveOutputs = false;

[H0,meanT] = centeredMTSensitivity(0,cfg.TuningSigmaDeg,true);
assert(abs(meanT-0.2611) < 0.002, 'Unexpected circular tuning mean.');
assert(H0 > 0, 'Centered drift sensitivity must be positive.');

bank = makeSparseTrialBank(6,100,cfg,cfg.Seed);
assert(all(abs(diag(bank.W)-H0) < 1e-12), ...
    'Every stream must be exactly aligned with one mechanism.');
assert(all(abs(bank.directionsDeg-(0:5)*60) < 1e-12), ...
    'Mechanisms are not evenly spaced or drift aligned.');
assert(max(bank.pedestalResponse)-min(bank.pedestalResponse) < 1e-12, ...
    'A common pedestal should add the same response to every mechanism.');

cfgPedestal = cfg;
cfgPedestal.StreamPedestal = 7.071;
bankPedestal = makeSparseTrialBank(6,100,cfgPedestal,cfg.Seed);
expectedShift = bankPedestal.pedestalResponse(1);
assert(max(abs((bankPedestal.responseChangedBase-bank.responseChangedBase) ...
    - expectedShift),[],'all') < 1e-10, ...
    'Uniform-pedestal invariance failed.');

results = runMaxPoolingNSweep(cfg);
matched = results.operatingTable(results.operatingTable.Mode=="Matched",:);
assert(all(abs(matched.Performance-cfg.TargetPerformance) < 0.04), ...
    'Matched-performance analysis trials are unexpectedly far from target.');
drift = results.summaryTable(results.summaryTable.OffsetDeg==0,:);
assert(all(abs(drift.NormalizedGain-1) < 1e-10), ...
    'Drift normalized gain must equal one.');
assert(all(isfinite(results.rawTable.Gain)), 'Nonfinite probit gain detected.');
assert(all(results.operatingTable.ProbitConverged), ...
    'At least one probit regression did not converge.');

close(results.figureHandles(ishandle(results.figureHandles)));
fprintf('validateSparseMaxPooling: all checks passed.\n');
end
