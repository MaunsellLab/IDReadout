function results = runPairedProbeTemporalMaxSweep(varargin)
% RUNPAIREDPROBETEMPORALMAXSWEEP  Difference of separate epoch maxima.
%
% HISTORICAL COMPARISON - NOT THE CURRENT MODEL
% ---------------------------------------------
% Winner-direction diagnostics showed that a 40% base drift virtually
% eliminates far-probe winners in this formulation and predicts a more
% negative 120-deg probe contribution at high base, opposite to IDR.  Use
% runPairedProbeDriftTemplateSweep for the current model.  This function is
% retained so the rejected computation remains reproducible.
%
% COMPLETE MODEL SPECIFICATION FOR THE TEMPORAL-MAX STAGE
% =======================================================
% Scientific question
% -------------------
% This stage asks whether change detection can preserve absolute-response
% competition.  Unlike the mechanism-wise delta model, direction pooling is
% performed separately in the prestep and step epochs and their scalar maxima
% are then subtracted.  For each spatial patch,
%
%   D_patch = max_m R_step,m - max_m R_pre,m,
%
% and the model chooses the changed patch when D_changed > D_noChange.  Thus
% the complete choice variable is a temporal difference within each patch
% followed by a spatial difference between patches.  Max pooling and temporal
% subtraction do not commute: this model is deliberately not
% max_m(R_step,m-R_pre,m).
%
% Stimulus and direction mechanisms
% ---------------------------------
% A dense circular bank of preferred directions uses 1-deg spacing by default.
% Drift is at 0 deg and paired probes are at +/-offset.  Each stream affects
% every mechanism through a centered wrapped-Gaussian MT sensitivity with
% sigma 37.5 deg.  Prestep absolute coherences are BaseDriftCoherence plus any
% DriftPedestal at 0 deg and ProbePedestal (7.071% by default) at both probes.
% The step retains these absolute coherences, adds one scalar Gaussian
% epoch-mean fluctuation per stream, and adds a fixed drift increment to the
% changed patch only.  Probe fluctuations are independent or yoked.  Changed
% and no-change patches are otherwise independent.
%
% Because the maxima operate on absolute landscapes, drift and probe pedestals
% determine competitive gaps and winner identities.  If the same mechanism
% wins both epochs, its constant baseline cancels from D_patch.  If the winner
% changes, pedestal geometry remains in the temporal difference.  This is the
% intended route to stable increment thresholds with preserved convergence and
% far-probe competition.
%
% Temporal sampling and internal noise
% ------------------------------------
% The first pass contains one integrated prestep sample and one integrated
% step sample; it does not simulate video frames or repeated working-memory
% updates.  Four independent smooth circular Gaussian fields represent the two
% epochs in the two spatial patches.  Their directional covariance is
%
%   rho(d) = exp{kappa[cos(d)-1]},
%
% with a 30-deg nominal correlation length and trial-wise circular-mean
% removal.  The dense grid therefore samples a correlated population field,
% not 360 independent competitors.  InternalStepNoiseSD is defined for the
% 250-ms step.  Assuming temporally white integrated noise, the 750-ms prestep
% SD is InternalStepNoiseSD*sqrt(T_step/T_pre).  Noise is homoscedastic across
% directions; the preceding Poisson stage showed that response-dependent count
% variance did not materially alter the gain curves.  There is no downstream
% noise or lapse process.
%
% One-time calibration and predictions
% ------------------------------------
% The step-epoch internal-noise SD is calibrated once so a 20% increment gives
% 75% correct at the low-base (7.071%), 90-deg, independent reference
% condition.  The same SD and duration scaling are then held fixed across all
% offsets, yoking modes, and base drifts (7.071% and 40% by default).  Every
% condition reports performance at the common 20% increment and a separately
% measured 75% threshold.  Stable thresholds outside the one calibration
% condition are predictions, not imposed operating points.
%
% Gain, winner, and transition diagnostics
% ----------------------------------------
% Correct/error choices are fit with simultaneous probit regression on the
% changed-minus-no-change scalar step fluctuations.  Independent probes have
% separate coefficients summarized by their mean and sum.  A yoked pair has
% one identifiable total coefficient; half is its per-stream equivalent.
%
% Prestep and step winners are classified as drift-associated,
% probe-associated, or other/intermediate.  The full 3-by-3 transition matrix
% is saved separately for changed and no-change patches.  Counterfactual
% perturbations further divide each stream's evidence derivative into:
% (1) trials whose temporal winners are identical, (2) trials whose prestep
% and step winners differ, and (3) additional step-winner switches caused by
% the perturbation itself.  These are mechanistic diagnostics, not fitted
% regression coefficients.
%
% A separate winner-direction diagnostic retains the full angular structure
% that is lost in the drift/probe/other categories.  Step winners are binned
% by preferred direction for the changed and no-change patches.  For each bin,
% the simulation saves (a) winner probability, (b) the centered-MT local
% sensitivity to the drift and to a per-stream paired-probe predictor, and
% (c) their product.  Summing the product over direction gives the ordinary
% winner-distribution-weighted fixed-winner derivative.  It therefore shows
% exactly which winner directions make a stream's contribution positive or
% negative.  It is not itself the fitted probit gain: perturbation-induced
% winner switches and unequal weighting of trials near the behavioral decision
% boundary remain separate effects in the counterfactual decomposition.
%
% Usage
% -----
%   results = runPairedProbeTemporalMaxSweep;
%   cfg = temporalMaxDefaults;
%   cfg.BaseDriftCoherences = [7.071 40];
%   results = runPairedProbeTemporalMaxSweep(cfg);

cfg=temporalMaxDefaults();
cfg=parseConfig(cfg,varargin{:}); pcfg=cfg.PairedConfig;
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end
calibrationBank=makeTemporalMaxTrialBank(pcfg.NCalibrationTrials,cfg, ...
    pcfg.Seed+410000);
analysisBank=makeTemporalMaxTrialBank(pcfg.NTrials,cfg,pcfg.Seed+420000);
if cfg.CalibrateInternalNoise
    calibration=calibrateTemporalMaxNoise(calibrationBank,cfg);
    cfg.InternalStepNoiseSD=calibration.internalStepNoiseSD;
else
    calibration=struct('internalStepNoiseSD',cfg.InternalStepNoiseSD, ...
        'performance',NaN,'performanceWithNoInternalNoise',NaN, ...
        'internalPreStepNoiseSD',cfg.InternalStepNoiseSD*sqrt( ...
        cfg.StepDurationSec/cfg.PreStepDurationSec), ...
        'baseDriftCoherence',cfg.ReferenceBaseDriftCoherence, ...
        'offsetDeg',pcfg.NoiseCalibrationOffsetDeg,'yoked', ...
        logical(pcfg.NoiseCalibrationYoked),'stepCoherence',pcfg.FixedStepCoherence);
end

nConditions=numel(cfg.BaseDriftCoherences)*numel(pcfg.OffsetsDeg)*2;
conditionPieces=cell(nConditions,1);
decompPieces=cell(2*nConditions,1);
transitionPieces=cell(18*nConditions,1);
winnerDirectionPieces=cell(nConditions,1);
piece=0; dd=0; tt=0;
for bb=1:numel(cfg.BaseDriftCoherences)
    base=cfg.BaseDriftCoherences(bb);
    for oo=1:numel(pcfg.OffsetsDeg)
        offset=pcfg.OffsetsDeg(oo);
        for yoked=[false true]
            piece=piece+1;
            threshold=calibrateTemporalMaxStep(calibrationBank,offset,yoked, ...
                base,cfg.InternalStepNoiseSD,cfg);
            sim=simulateTemporalMaxCondition(analysisBank,offset,yoked, ...
                pcfg.FixedStepCoherence,base,cfg.InternalStepNoiseSD,cfg);
            conditionPieces{piece}=conditionTable(sim,threshold,cfg);
            winnerDirectionPieces{piece}=winnerDirectionTable(sim,cfg);
            for streamNumber=1:2
                dd=dd+1; decompPieces{dd}=decompositionTable(sim,streamNumber,cfg);
            end
            for patchNumber=1:2
                matrix=sim.transitions.changed;
                patchName="Changed";
                if patchNumber==2
                    matrix=sim.transitions.noChange; patchName="NoChange";
                end
                for from=1:3
                    for to=1:3
                        tt=tt+1;
                        transitionPieces{tt}=transitionTable(sim,patchName, ...
                            from,to,matrix(from,to),cfg);
                    end
                end
            end
        end
    end
end

results=struct(); results.config=cfg; results.noiseCalibration=calibration;
results.noiseCalibrationTable=struct2table(calibration);
results.conditionTable=vertcat(conditionPieces{:});
results.decompositionTable=vertcat(decompPieces{:});
results.transitionTable=vertcat(transitionPieces{:});
results.winnerDirectionTable=vertcat(winnerDirectionPieces{:});
if cfg.CreateFigures
    results.figureHandles=plotPairedProbeTemporalMaxSweep(results);
else
    results.figureHandles=gobjects(0);
end
if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_TemporalMax_OffsetSummary.csv'));
    writetable(results.decompositionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_TemporalMax_Decomposition.csv'));
    writetable(results.transitionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_TemporalMax_WinnerTransitions.csv'));
    writetable(results.winnerDirectionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_TemporalMax_WinnerDirectionContributions.csv'));
    writetable(results.noiseCalibrationTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_TemporalMax_NoiseCalibration.csv'));
    handles=results.figureHandles; results.figureHandles=[];
    save(fullfile(cfg.OutputDir,'PairedProbe_TemporalMax_Results.mat'), ...
        'results','-v7.3'); results.figureHandles=handles;
    if cfg.CreateFigures
        savePairedProbeTemporalMaxFigures(handles,cfg.OutputDir);
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
validateattributes(cfg.BaseDriftCoherences,{'numeric'},{'vector','nonnegative'});
validateattributes(cfg.PreStepDurationSec,{'numeric'},{'scalar','positive'});
validateattributes(cfg.StepDurationSec,{'numeric'},{'scalar','positive'});
validateattributes(cfg.WinnerDirectionBinWidthDeg,{'numeric'}, ...
    {'scalar','positive','<=',360});
if abs(360/cfg.WinnerDirectionBinWidthDeg- ...
        round(360/cfg.WinnerDirectionBinWidthDeg))>1e-10
    error('WinnerDirectionBinWidthDeg must divide 360 exactly.');
end
validateattributes(cfg.WinnerDirectionPlotOffsetDeg,{'numeric'}, ...
    {'scalar','nonnegative','<=',180});
end

function T=conditionTable(sim,threshold,cfg)
g=sim.gain; w=sim.winners; tr=sim.transitions; pcfg=cfg.PairedConfig;
T=table(sim.baseDriftCoherence,sim.offsetDeg,string(sim.yoking), ...
    sim.performance,sim.stepCoherence,threshold,g.drift,g.driftSE, ...
    g.probePerStream,g.probePerStreamSE,g.probeTotal,g.probeTotalSE, ...
    g.probePerStream/g.drift,w.changedDrift,w.changedProbe,w.changedOther, ...
    w.noChangeDrift,w.noChangeProbe,w.noChangeOther,tr.changedSame,tr.noChangeSame, ...
    sim.fit.intercept,sim.fit.logLikelihood,sim.fit.converged, ...
    cfg.InternalStepNoiseSD,cfg.InternalStepNoiseSD*sqrt( ...
    cfg.StepDurationSec/cfg.PreStepDurationSec),pcfg.ProbePedestal, ...
    pcfg.InternalCorrelationLengthDeg,circularNoiseEffectiveRank( ...
    round(360/pcfg.BankSpacingDeg),pcfg), ...
    'VariableNames',{'BaseDriftCoherence','OffsetDeg','Yoking','Performance', ...
    'StepCoherence','Threshold75','DriftGain','DriftGainSE', ...
    'ProbeGainPerStream','ProbeGainPerStreamSE','ProbeGainTotal', ...
    'ProbeGainTotalSE','NormalizedProbePerStream','ChangedDriftWinner', ...
    'ChangedProbeWinner','ChangedOtherWinner','NoChangeDriftWinner', ...
    'NoChangeProbeWinner','NoChangeOtherWinner','ChangedTemporalSameWinner', ...
    'NoChangeTemporalSameWinner','ProbitIntercept','LogLikelihood', ...
    'ProbitConverged','InternalStepNoiseSD','InternalPreStepNoiseSD', ...
    'ProbePedestal','InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end

function T=decompositionTable(sim,streamNumber,cfg)
if streamNumber==1, name="Drift"; d=sim.decomposition.drift;
else, name="ProbePerStream"; d=sim.decomposition.probePerStream; end
T=table(sim.baseDriftCoherence,sim.offsetDeg,string(sim.yoking),name, ...
    d.decisionTemporalSameFixed,d.decisionTemporalDifferentFixed, ...
    d.decisionPerturbationSwitch,d.decisionTotal, ...
    d.changedTemporalSameProbability,d.noChangeTemporalSameProbability, ...
    d.changedPerturbationSwitchProbability,d.noChangePerturbationSwitchProbability, ...
    cfg.InternalStepNoiseSD,'VariableNames',{'BaseDriftCoherence','OffsetDeg', ...
    'Yoking','Stream','TemporalSameWinner','TemporalDifferentWinner', ...
    'PerturbationWinnerSwitch','DecisionTotal','ChangedTemporalSameProbability', ...
    'NoChangeTemporalSameProbability','ChangedPerturbationSwitchProbability', ...
    'NoChangePerturbationSwitchProbability','InternalStepNoiseSD'});
end

function T=transitionTable(sim,patchName,from,to,fraction,cfg)
names=["Drift","Probe","Other"];
T=table(sim.baseDriftCoherence,sim.offsetDeg,string(sim.yoking),patchName, ...
    names(from),names(to),fraction,cfg.InternalStepNoiseSD, ...
    'VariableNames',{'BaseDriftCoherence','OffsetDeg','Yoking','Patch', ...
    'PreStepWinner','StepWinner','Fraction','InternalStepNoiseSD'});
end

function T=winnerDirectionTable(sim,cfg)
d=sim.winnerDirection; n=numel(d.directionDeg);
T=table(repmat(sim.baseDriftCoherence,n,1),repmat(sim.offsetDeg,n,1), ...
    repmat(string(sim.yoking),n,1),d.directionDeg, ...
    d.changedWinnerFraction,d.noChangeWinnerFraction, ...
    d.patchAverageWinnerFraction,d.driftLocalSensitivity, ...
    d.probeLocalSensitivityPerStream,d.driftWinnerWeightedContribution, ...
    d.probeWinnerWeightedContributionPerStream, ...
    repmat(d.driftFixedWinnerTotal,n,1), ...
    repmat(d.probeFixedWinnerTotalPerStream,n,1), ...
    repmat(cfg.WinnerDirectionBinWidthDeg,n,1), ...
    'VariableNames',{'BaseDriftCoherence','OffsetDeg','Yoking', ...
    'WinnerDirectionDeg','ChangedWinnerFraction','NoChangeWinnerFraction', ...
    'PatchAverageWinnerFraction','DriftLocalSensitivity', ...
    'ProbeLocalSensitivityPerStream','DriftWinnerWeightedContribution', ...
    'ProbeWinnerWeightedContributionPerStream','DriftFixedWinnerTotal', ...
    'ProbeFixedWinnerTotalPerStream','DirectionBinWidthDeg'});
end
