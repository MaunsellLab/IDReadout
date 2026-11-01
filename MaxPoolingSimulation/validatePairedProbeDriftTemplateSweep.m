function validatePairedProbeDriftTemplateSweep()
% VALIDATEPAIREDPROBEDRIFTTEMPLATESWEEP  Fast invariance and output checks.

cfg=driftTemplateDefaults();
cfg.PairedConfig.OffsetsDeg=[120 180];
cfg.CancellationFractions=[0 0.85 1];
cfg.PairedConfig.BankSpacingDeg=2;
cfg.WinnerDirectionBinWidthDeg=10;
cfg.PairedConfig.NTrials=2000;
cfg.PairedConfig.NCalibrationTrials=3000;
cfg.PairedConfig.CalibrationTolerance=0.015;
cfg.PairedConfig.CalibrationMaxIterations=22;
cfg.CreateFigures=false; cfg.SaveOutputs=false; cfg.Visible='off';
results=runPairedProbeDriftTemplateSweep(cfg);
T=results.conditionTable; W=results.winnerDirectionTable;
nConditions=2*2*2*3;
assert(height(T)==nConditions,'Unexpected number of conditions.');
assert(height(W)==nConditions*(360/cfg.WinnerDirectionBinWidthDeg), ...
    'Unexpected number of winner-direction rows.');
assert(all(isfinite(T.Threshold75)) && all(T.Threshold75>0), ...
    'Invalid threshold.');
assert(all(T.ProbitConverged),'At least one probit fit did not converge.');
groups=findgroups(W.BaseDriftCoherence,W.CancellationFraction, ...
    W.OffsetDeg,W.Yoking);
sums=splitapply(@sum,W.PatchAverageWinnerFraction,groups);
assert(max(abs(sums-1))<1e-10, ...
    'Winner-direction distributions do not sum to one.');

full=T(T.CancellationFraction==1,:);
for mode=["Independent","Yoked"]
    for offset=[120 180]
        S=sortrows(full(full.Yoking==mode & full.OffsetDeg==offset,:), ...
            'BaseDriftCoherence');
        vars={'Performance','Threshold75','DriftGain','ProbeGainPerStream', ...
            'ChangedProbeWinner','NoChangeProbeWinner'};
        for k=1:numel(vars)
            assert(abs(S.(vars{k})(1)-S.(vars{k})(2))<1e-10, ...
                'Full cancellation should remove base-drift dependence.');
        end
    end
end

S=T(T.OffsetDeg==120 & T.Yoking=="Independent" & ...
    T.CancellationFraction==0,:);
S=sortrows(S,'BaseDriftCoherence');
assert(S.ChangedProbeWinner(2)<S.ChangedProbeWinner(1), ...
    'The uncancelled high base should suppress probe winners.');
fprintf('Drift-template validation passed. Internal-noise SD = %.4g.\n', ...
    results.noiseCalibration.internalNoiseSD);
end
