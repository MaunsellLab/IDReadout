function validateSplitSignalModelRecovery()
% VALIDATESPLITSIGNALMODELRECOVERY  Small end-to-end recovery test.

cfg=splitSignalRecoveryDefaults();
cfg.NCalibrationTrials=800; cfg.NStimulusTrials=450;
cfg.NLinkCalibrationTrials=450;
cfg.NRecoveryReplicates=2; cfg.NGainReplicates=2;
cfg.TrialsPerSession=100; cfg.MaxRecoverySessions=2;
cfg.ThresholdMaxIterations=5; cfg.ThresholdTolerance=2;
cfg.CreateFigures=false; cfg.SaveOutputs=false;
cfg.Verbose=false;
pairs=structToPairs(cfg); results=runSplitSignalModelRecovery(pairs{:});
assert(all(isfinite(results.thresholds(:))));
assert(all(size(results.thresholds)==[3 3]));
for m=1:3
    assert(all(size(results.recovery(m).probability)==[2 3]));
    assert(all(abs(sum(results.recovery(m).probability,2)-1)<1e-9));
    assert(numel(results.gains(m).stream1Mean)==3);
end
fprintf('validateSplitSignalModelRecovery passed.\n');
end

function pairs=structToPairs(s)
names=fieldnames(s); pairs=cell(1,2*numel(names));
for k=1:numel(names), pairs{2*k-1}=names{k}; pairs{2*k}=s.(names{k}); end
end
