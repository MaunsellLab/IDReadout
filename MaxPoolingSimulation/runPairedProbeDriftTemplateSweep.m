function results = runPairedProbeDriftTemplateSweep(varargin)
% RUNPAIREDPROBEDRIFTTEMPLATESWEEP  Cancel perceived drift before max pooling.
% HISTORICAL INTERMEDIATE STAGE: the current biologically constrained model is
% runPairedProbeAdaptationSweep.  This runner is retained to show what exact or
% partial drift-template cancellation does without firing-dependent adaptation.
%
% COMPLETE MODEL SPECIFICATION FOR THE DRIFT-TEMPLATE STAGE
% =========================================================
% Scientific purpose
% ------------------
% The preceding absolute-response temporal-max model applied the direction
% maximum separately to prestep and step responses and then subtracted the two
% maxima.  At a 40% base drift, that model made the drift mechanism win almost
% every step-period competition, eliminating 120-deg probe winners and making
% the probe contribution strongly negative.  This is opposite to the positive
% far-probe gains observed in high-base IDR.  The absolute-response temporal
% max is therefore not part of the present model; its code remains in the
% package only as a historical comparison.
%
% Perceptual drift-template cancellation
% --------------------------------------
% A dense circular bank of direction mechanisms receives the absolute step
% response after subtraction of an additive template for the perceived
% prestep drift.  For mechanism m in either spatial patch,
%
%   E_m = sum_j C_step,j H(theta_j-phi_m)
%         - alpha B_drift H(0-phi_m) + epsilon_m,
%
% where B_drift is BaseDriftCoherence+DriftPedestal and alpha is the explicit
% CancellationFraction.  Equivalently, the residual sustained drift entering
% the max is (1-alpha)B_drift.  Alpha=0 leaves the absolute step landscape;
% alpha=1 removes the perceived drift template completely.  The model then
% computes max_m E_m separately in the changed and no-change patches and
% chooses the changed patch when its maximum is larger.  There is no prestep
% maximum and no subtraction of scalar epoch maxima.
%
% The cancellation is selective for the perceived drift template.  The two
% 7.071% probe pedestals are intentionally retained, as are the scalar drift
% and probe coherence fluctuations and the fixed positive drift increment.
% This distinguishes the model from pure mechanism-wise delta coherence, in
% which all constant prestep components cancel.  The operation is additive,
% not divisive, so it does not impose Weber-fraction scaling.
%
% Stimulus streams and direction bank
% -----------------------------------
% Drift is at 0 deg and probes are at +/-offset.  Preferred directions span
% the full circle at 1-deg spacing by default.  Each stream affects every
% mechanism through centered wrapped-Gaussian MT sensitivity with sigma 37.5
% deg.  The changed patch alone receives a fixed 20% drift increment.  One
% scalar Gaussian epoch-mean fluctuation is drawn per stream and patch.  Probe
% fluctuations are independent or yoked; changed and no-change patches remain
% independent.  The abstraction does not simulate video frames or enforce a
% fixed dot budget.
%
% Internal noise and calibration
% ------------------------------
% Additive mechanism noise is a smooth circular Gaussian field with the same
% 30-deg nominal correlation length and trial-wise circular-mean removal used
% in the preceding dense-bank stages.  It is not 360 independent samples.
% There is no downstream noise or lapse process.
%
% The internal-noise SD is calibrated once so the 20% increment produces 75%
% correct at the low-base, 90-deg, independent-probe reference condition with
% ReferenceCancellationFraction (1 by default).  That SD is held fixed across
% every base drift, offset, yoking condition, and cancellation fraction.  Each
% condition reports performance at the common 20% increment and a separately
% measured 75% threshold.  These are predictions outside the one calibration
% condition.
%
% Cancellation sweep and gain recovery
% -------------------------------------
% The default focused sweep uses offsets 120, 170, and 180 deg, base drifts
% 7.071% and 40%, both yoking modes, and alpha values
% [0 .5 .75 .825 .875 .925 .975 1].  It is designed to locate a range that
% simultaneously restores high-base far-probe winners, retains low-base
% negative 120-deg gain, preserves the 170-180-deg rebound, and keeps the 20%
% threshold approximately stable.  A later full-offset sweep can use the same
% runner after editing PairedConfig.OffsetsDeg.
%
% Correct/error choices are fit with simultaneous probit regression on the
% changed-minus-no-change stream fluctuations.  Independent probes receive
% separate coefficients summarized by their mean and sum.  A yoked pair has
% one identifiable total coefficient; half is reported as its per-stream
% equivalent.  Counterfactual perturbations divide evidence derivatives into
% fixed-winner and perturbation-induced winner-switch components.
%
% Winner-direction diagnostics
% ----------------------------
% For every condition, actual changed- and no-change-patch winners are binned
% by preferred direction.  The output saves winner probabilities, local drift
% and per-stream probe sensitivities, and their signed products.  Summing the
% products gives the ordinary winner-distribution-weighted fixed-winner
% derivative.  This diagnostic does not include perturbation-induced switches
% or extra weighting of trials near the behavioral decision boundary and is
% therefore not identical to the recovered probit gain.
%
% Usage
% -----
%   results = runPairedProbeDriftTemplateSweep;
%   cfg = driftTemplateDefaults;
%   cfg.CancellationFractions = [0.75 0.85 0.95 1];
%   results = runPairedProbeDriftTemplateSweep(cfg);

cfg=driftTemplateDefaults();
cfg=parseConfig(cfg,varargin{:}); pcfg=cfg.PairedConfig;
pcfg.EffectiveNoiseRank=circularNoiseEffectiveRank( ...
    round(360/pcfg.BankSpacingDeg),pcfg);
cfg.PairedConfig=pcfg;
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

calibrationBank=makeDensePairedTrialBank(pcfg.NCalibrationTrials,pcfg, ...
    pcfg.Seed+510000);
analysisBank=makeDensePairedTrialBank(pcfg.NTrials,pcfg,pcfg.Seed+520000);
if pcfg.CalibrateInternalNoise
    calibration=calibrateDriftTemplateNoise(calibrationBank,cfg);
    pcfg.InternalNoiseSD=calibration.internalNoiseSD;
    cfg.PairedConfig=pcfg;
else
    calibration=struct('internalNoiseSD',pcfg.InternalNoiseSD, ...
        'performance',NaN,'performanceWithNoInternalNoise',NaN, ...
        'baseDriftCoherence',cfg.ReferenceBaseDriftCoherence, ...
        'cancellationFraction',cfg.ReferenceCancellationFraction, ...
        'offsetDeg',pcfg.NoiseCalibrationOffsetDeg, ...
        'yoked',logical(pcfg.NoiseCalibrationYoked), ...
        'stepCoherence',pcfg.FixedStepCoherence);
end
calibration.correlationLengthDeg=pcfg.InternalCorrelationLengthDeg;
calibration.effectiveNoiseRank=pcfg.EffectiveNoiseRank;
calibration.removeCircularMean=logical(pcfg.RemoveInternalCircularMean);

nConditions=numel(cfg.BaseDriftCoherences)*numel(pcfg.OffsetsDeg)*2* ...
    numel(cfg.CancellationFractions);
conditionPieces=cell(nConditions,1);
decompositionPieces=cell(2*nConditions,1);
winnerDirectionPieces=cell(nConditions,1);
piece=0; dd=0;
for bb=1:numel(cfg.BaseDriftCoherences)
    base=cfg.BaseDriftCoherences(bb);
    for oo=1:numel(pcfg.OffsetsDeg)
        offset=pcfg.OffsetsDeg(oo);
        for yoked=[false true]
            for aa=1:numel(cfg.CancellationFractions)
                alpha=cfg.CancellationFractions(aa); piece=piece+1;
                threshold=calibrateDriftTemplateStep(calibrationBank,offset, ...
                    yoked,base,alpha,pcfg.InternalNoiseSD,cfg);
                sim=simulateDriftTemplateCondition(analysisBank,offset,yoked, ...
                    pcfg.FixedStepCoherence,base,alpha,pcfg.InternalNoiseSD,cfg);
                conditionPieces{piece}=conditionTable(sim,threshold,cfg);
                winnerDirectionPieces{piece}=winnerDirectionTable(sim,cfg);
                for streamNumber=1:2
                    dd=dd+1;
                    decompositionPieces{dd}=decompositionTable(sim, ...
                        streamNumber,cfg);
                end
            end
        end
    end
end

results=struct(); results.config=cfg; results.noiseCalibration=calibration;
results.noiseCalibrationTable=struct2table(calibration);
results.conditionTable=vertcat(conditionPieces{:});
results.decompositionTable=vertcat(decompositionPieces{:});
results.winnerDirectionTable=vertcat(winnerDirectionPieces{:});
if cfg.CreateFigures
    results.figureHandles=plotPairedProbeDriftTemplateSweep(results);
else
    results.figureHandles=gobjects(0);
end
if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_DriftTemplate_ConditionSummary.csv'));
    writetable(results.decompositionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_DriftTemplate_Decomposition.csv'));
    writetable(results.winnerDirectionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_DriftTemplate_WinnerDirectionContributions.csv'));
    writetable(results.noiseCalibrationTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_DriftTemplate_NoiseCalibration.csv'));
    handles=results.figureHandles; results.figureHandles=[];
    save(fullfile(cfg.OutputDir,'PairedProbe_DriftTemplate_Results.mat'), ...
        'results','-v7.3'); results.figureHandles=handles;
    if cfg.CreateFigures
        savePairedProbeDriftTemplateFigures(handles,cfg.OutputDir);
    end
end
end

function cfg=parseConfig(cfg,varargin)
if isempty(varargin)
elseif numel(varargin)==1 && isstruct(varargin{1})
    supplied=varargin{1}; names=fieldnames(supplied);
    for k=1:numel(names)
        if ~isfield(cfg,names{k}), error('Unknown configuration field: %s',names{k}); end
        cfg.(names{k})=supplied.(names{k});
    end
elseif mod(numel(varargin),2)==0
    for k=1:2:numel(varargin)
        name=char(varargin{k});
        if ~isfield(cfg,name), error('Unknown configuration field: %s',name); end
        cfg.(name)=varargin{k+1};
    end
else
    error('Use a configuration struct or name/value pairs.');
end
pcfg=cfg.PairedConfig;
validateattributes(cfg.BaseDriftCoherences,{'numeric'}, ...
    {'vector','nonnegative'});
validateattributes(cfg.CancellationFractions,{'numeric'}, ...
    {'vector','>=',0,'<=',1});
validateattributes(cfg.ReferenceBaseDriftCoherence,{'numeric'}, ...
    {'scalar','nonnegative'});
validateattributes(cfg.ReferenceCancellationFraction,{'numeric'}, ...
    {'scalar','>=',0,'<=',1});
validateattributes(cfg.WinnerDirectionBinWidthDeg,{'numeric'}, ...
    {'scalar','positive','<=',360});
if abs(360/cfg.WinnerDirectionBinWidthDeg- ...
        round(360/cfg.WinnerDirectionBinWidthDeg))>1e-10
    error('WinnerDirectionBinWidthDeg must divide 360 exactly.');
end
validateattributes(cfg.WinnerDirectionPlotOffsetDeg,{'numeric'}, ...
    {'scalar','>=',0,'<=',180});
validateattributes(cfg.WinnerDirectionPlotBase,{'numeric'}, ...
    {'scalar','nonnegative'});
validateattributes(cfg.WinnerDirectionPlotFractions,{'numeric'}, ...
    {'vector','>=',0,'<=',1});
validateattributes(pcfg.OffsetsDeg,{'numeric'},{'vector','>=',0,'<=',180});
validateattributes(pcfg.BankSpacingDeg,{'numeric'}, ...
    {'scalar','positive','<=',10});
if abs(360/pcfg.BankSpacingDeg-round(360/pcfg.BankSpacingDeg))>1e-10
    error('BankSpacingDeg must divide 360 exactly.');
end
end

function T=conditionTable(sim,threshold,cfg)
g=sim.gain; w=sim.winners; pcfg=cfg.PairedConfig;
T=table(sim.baseDriftCoherence,sim.cancellationFraction, ...
    sim.residualDriftCoherence,sim.cancelledDriftCoherence,sim.offsetDeg, ...
    string(sim.yoking),sim.performance,sim.stepCoherence,threshold, ...
    g.drift,g.driftSE,g.probePerStream,g.probePerStreamSE,g.probeTotal, ...
    g.probeTotalSE,g.probePerStream/g.drift,w.changedDrift,w.changedProbe, ...
    w.changedOther,w.noChangeDrift,w.noChangeProbe,w.noChangeOther, ...
    sim.fit.intercept,sim.fit.logLikelihood,sim.fit.converged, ...
    pcfg.InternalNoiseSD,pcfg.ProbePedestal, ...
    pcfg.InternalCorrelationLengthDeg,pcfg.EffectiveNoiseRank, ...
    'VariableNames',{'BaseDriftCoherence','CancellationFraction', ...
    'ResidualDriftCoherence','CancelledDriftCoherence','OffsetDeg','Yoking', ...
    'Performance','StepCoherence','Threshold75','DriftGain','DriftGainSE', ...
    'ProbeGainPerStream','ProbeGainPerStreamSE','ProbeGainTotal', ...
    'ProbeGainTotalSE','NormalizedProbePerStream','ChangedDriftWinner', ...
    'ChangedProbeWinner','ChangedOtherWinner','NoChangeDriftWinner', ...
    'NoChangeProbeWinner','NoChangeOtherWinner','ProbitIntercept', ...
    'LogLikelihood','ProbitConverged','InternalNoiseSD','ProbePedestal', ...
    'InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end

function T=decompositionTable(sim,streamNumber,cfg)
if streamNumber==1, name="Drift"; d=sim.decomposition.drift;
else, name="ProbePerStream"; d=sim.decomposition.probePerStream; end
T=table(sim.baseDriftCoherence,sim.cancellationFraction,sim.offsetDeg, ...
    string(sim.yoking),name,d.changedFixed,d.changedSwitch,d.noChangeFixed, ...
    d.noChangeSwitch,d.decisionFixed,d.decisionSwitch,d.decisionTotal, ...
    d.changedSwitchProbability,d.noChangeSwitchProbability, ...
    cfg.PairedConfig.InternalNoiseSD, ...
    'VariableNames',{'BaseDriftCoherence','CancellationFraction','OffsetDeg', ...
    'Yoking','Stream','ChangedFixedWinner','ChangedWinnerSwitch', ...
    'NoChangeFixedWinner','NoChangeWinnerSwitch','DecisionFixedWinner', ...
    'DecisionWinnerSwitch','DecisionTotal','ChangedSwitchProbability', ...
    'NoChangeSwitchProbability','InternalNoiseSD'});
end

function T=winnerDirectionTable(sim,cfg)
d=sim.winnerDirection; n=numel(d.directionDeg);
T=table(repmat(sim.baseDriftCoherence,n,1), ...
    repmat(sim.cancellationFraction,n,1),repmat(sim.offsetDeg,n,1), ...
    repmat(string(sim.yoking),n,1),d.directionDeg, ...
    d.changedWinnerFraction,d.noChangeWinnerFraction, ...
    d.patchAverageWinnerFraction,d.driftLocalSensitivity, ...
    d.probeLocalSensitivityPerStream,d.driftWinnerWeightedContribution, ...
    d.probeWinnerWeightedContributionPerStream, ...
    repmat(d.driftFixedWinnerTotal,n,1), ...
    repmat(d.probeFixedWinnerTotalPerStream,n,1), ...
    repmat(cfg.WinnerDirectionBinWidthDeg,n,1), ...
    'VariableNames',{'BaseDriftCoherence','CancellationFraction','OffsetDeg', ...
    'Yoking','WinnerDirectionDeg','ChangedWinnerFraction', ...
    'NoChangeWinnerFraction','PatchAverageWinnerFraction', ...
    'DriftLocalSensitivity','ProbeLocalSensitivityPerStream', ...
    'DriftWinnerWeightedContribution', ...
    'ProbeWinnerWeightedContributionPerStream','DriftFixedWinnerTotal', ...
    'ProbeFixedWinnerTotalPerStream','DirectionBinWidthDeg'});
end
