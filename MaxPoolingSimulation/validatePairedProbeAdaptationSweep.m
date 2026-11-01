function validatePairedProbeAdaptationSweep()
% VALIDATEPAIREDPROBEADAPTATIONSWEEP  Fast structural and invariance checks.
cfg=adaptationSweepDefaults();
cfg.PairedConfig.OffsetsDeg=[120 170]; cfg.AftereffectAt60Values=[0 4 7];
cfg.PairedConfig.BankSpacingDeg=2; cfg.WinnerDirectionBinWidthDeg=10;
cfg.PairedConfig.NTrials=2000; cfg.PairedConfig.NCalibrationTrials=3000;
cfg.PairedConfig.CalibrationTolerance=0.015; cfg.PairedConfig.CalibrationMaxIterations=22;
cfg.CreateFigures=false; cfg.SaveOutputs=false; cfg.Visible='off';
r=runPairedProbeAdaptationSweep(cfg); T=r.conditionTable; W=r.winnerDirectionTable;
n=2*2*2*3; assert(height(T)==n,'Unexpected condition count.');
assert(height(r.baseContrastTable)==n/2,'Unexpected base-contrast count.');
assert(height(W)==n*(360/cfg.WinnerDirectionBinWidthDeg),'Unexpected winner-row count.');
assert(all(isfinite(T.Threshold75)) && all(T.Threshold75>0),'Invalid threshold.');
assert(all(T.ProbitConverged),'At least one probit fit did not converge.');
groups=findgroups(W.BaseDriftCoherence,W.AftereffectAt60,W.OffsetDeg,W.Yoking);
sums=splitapply(@sum,W.PatchAverageWinnerFraction,groups);
assert(max(abs(sums-1))<1e-10,'Winner-direction distributions do not sum to one.');

% With zero adaptation, exact drift cancellation makes the two bases identical.
Z=T(T.AftereffectAt60==0,:);
for mode=["Independent","Yoked"], for offset=[120 170]
    S=sortrows(Z(Z.Yoking==mode & Z.OffsetDeg==offset,:),'BaseDriftCoherence');
    vars={'Performance','Threshold75','DriftGain','ProbeGainPerStream', ...
        'ChangedProbeWinner','NoChangeProbeWinner'};
    for k=1:numel(vars)
        assert(abs(S.(vars{k})(1)-S.(vars{k})(2))<1e-10, ...
            'Zero adaptation should remove base-drift dependence.');
    end
end, end

reference=T(T.BaseDriftCoherence==cfg.ReferenceBaseDriftCoherence & ...
    T.OffsetDeg==cfg.PairedConfig.NoiseCalibrationOffsetDeg & ...
    T.Yoking=="Independent" & ...
    T.AftereffectAt60==cfg.ReferenceAftereffectAt60,:);
assert(height(reference)==1,'Reference calibration condition is not unique.');
assert(abs(reference.Threshold75-cfg.PairedConfig.FixedStepCoherence)<1.5, ...
    'Measured reference threshold is inconsistent with the 15.3% anchor.');

% The implemented field must scale with base and suppress the drift center.
S=sortrows(T(T.AftereffectAt60==7 & T.OffsetDeg==120 & ...
    T.Yoking=="Independent",:),'BaseDriftCoherence');
assert(S.ActualAftereffect(2)>S.ActualAftereffect(1), ...
    'Adaptation did not scale with prestep drift.');
assert(S.NoChangeDriftWinner(2)<S.NoChangeDriftWinner(1), ...
    'High-base adaptation should reduce no-change drift winners.');
contrast=r.baseContrastTable(r.baseContrastTable.OffsetDeg==120 & ...
    r.baseContrastTable.Yoking=="Independent",:);
assert(all(isfinite(contrast.ThresholdIncrease)), ...
    'Base-dependent threshold increments were not computed.');
fprintf('Adaptation validation passed. Internal-noise SD = %.4g.\n', ...
    r.noiseCalibration.internalNoiseSD);
end
