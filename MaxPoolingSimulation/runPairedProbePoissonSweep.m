function results = runPairedProbePoissonSweep(varargin)
% RUNPAIREDPROBEPOISSONSWEEP  Max pooling with response-dependent count noise.
%
% COMPLETE MODEL SPECIFICATION FOR THE POOLED-POISSON VARIANCE STAGE
% =================================================================
% Purpose and stimulus representation
% -----------------------------------
% This stage retains the paired-probe change detector and asks whether the
% variance expected from neuronal spike counts changes the max-pooling
% predictions.  The stimulus has drift at 0 deg and probe streams at
% +/-offset.  Each eligible direction mechanism first computes its difference
% between mean step-epoch and mean prestep-epoch activity.  There is no
% frame-resolved coherence or temporal filtering in this first pass.
%
% Direction bank and mean response
% --------------------------------
% Preferred directions form a dense circular bank (1 deg by default).  Stream
% coherence affects each mechanism through the same centered Gaussian MT
% sensitivity used in the preceding simulations, with sigma 37.5 deg.  The
% tuning scale is normalized so an aligned 1%-coherence increment changes the
% expected firing rate by AlignedRateSlopeSpPerSecPerPct (default 0.5 sp/s).
% Random motion contributes RandomMotionRateSpPerSec (default 40 sp/s).
%
% Absolute prestep rates include BaseDriftCoherence at 0 deg and the 7.071%
% pedestal at each probe.  Step rates retain those means, add the scalar
% Gaussian epoch-mean changes already used for the drift and probe streams,
% and add a fixed drift increment to the changed patch only.  Mean decision
% evidence is the mechanism-wise step-minus-prestep rate difference, expressed
% back in coherence-response units.  Thus constant bases cancel from the mean
% but still affect spike-count variance.
%
% Pooled-Poisson variance
% -----------------------
% Conditional count noise is approximated as Gaussian with Poisson/Fano
% variance.  For mechanism m,
%
%   Var(delta_m) = Fano/N_eff * [r_step,m/T_step + r_pre,m/T_pre] / a^2,
%
% where a is the rate-per-model-unit scale and N_eff is one fitted effective
% independent-neuron pool size.  Prestep and step fields are independent;
% changed and no-change patches are independent.  Within each field, direction
% samples have the same smooth circular 30-deg correlation used previously,
% with trial-wise circular-mean removal.  The dense 1-deg grid therefore
% samples a correlated population field rather than creating 360 independent
% competitors.  Rates are bounded below only to keep Poisson variance valid,
% and clamp fractions are reported.  There is no downstream decision noise or
% lapse process.
%
% Calibration and tests
% ---------------------
% N_eff is calibrated once so a fixed 20% drift increment gives 75% correct at
% the low-base (7.071%), 90-deg, independent-probe reference condition.  That
% same pool size is then held fixed across all offsets, yoking conditions, and
% base drifts (7.071% and 40% by default).  Each condition reports performance
% at the common 20% increment and separately measures the increment producing
% 75% correct.  The decisive test is whether these thresholds remain near 20%
% across base drift without erasing the max-competition gain effects.
%
% Gains and winners
% -----------------
% Correct/error choices are fit by simultaneous probit regression on the
% change-minus-no-change scalar stream fluctuations.  Independent probes have
% separate coefficients; their mean is the per-stream gain and their sum the
% pair gain.  A yoked pair has one identifiable coefficient, with half reported
% as its per-stream equivalent.  Winning mechanisms are classified as drift,
% probe, or other/intermediate.  Optional dashed overlays reproduce the prior
% homoscedastic Gaussian-noise model for comparison.
%
% Usage
% -----
%   results = runPairedProbePoissonSweep;
%   cfg = pairedProbePoissonDefaults;
%   cfg.BaseDriftCoherences = [7.071 20 40];
%   results = runPairedProbePoissonSweep(cfg);

cfg=pairedProbePoissonDefaults();
cfg=parseConfig(cfg,varargin{:});
pcfg=cfg.PairedConfig;
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

calibrationBank=makePoissonPairedTrialBank(pcfg.NCalibrationTrials,cfg, ...
    pcfg.Seed+310000);
analysisBank=makePoissonPairedTrialBank(pcfg.NTrials,cfg,pcfg.Seed+320000);
poolCalibration=calibratePoissonPoolSize(calibrationBank,cfg);
effectivePoolSize=poolCalibration.effectivePoolSize;

nConditions=numel(cfg.BaseDriftCoherences)*numel(pcfg.OffsetsDeg)*2;
pieces=cell(nConditions,1); piece=0;
for bb=1:numel(cfg.BaseDriftCoherences)
    base=cfg.BaseDriftCoherences(bb);
    for oo=1:numel(pcfg.OffsetsDeg)
        offset=pcfg.OffsetsDeg(oo);
        for yoked=[false true]
            piece=piece+1;
            threshold=calibratePoissonPairedProbeStep(calibrationBank, ...
                offset,yoked,base,effectivePoolSize,cfg);
            sim=simulatePoissonPairedProbeCondition(analysisBank,offset,yoked, ...
                pcfg.FixedStepCoherence,base,effectivePoolSize,cfg);
            pieces{piece}=makeConditionTable(sim,threshold,cfg);
        end
    end
end

homoscedasticTable=table();
if cfg.IncludeHomoscedasticReference
    referenceCfg=pcfg;
    referenceCfg.SaveOutputs=false;
    referenceCfg.CreateFigures=false;
    referenceCfg.Visible='off';
    reference=runPairedProbeOffsetSweep(referenceCfg);
    homoscedasticTable=reference.conditionTable;
end

results=struct();
results.config=cfg;
results.poolCalibration=poolCalibration;
results.poolCalibrationTable=struct2table(poolCalibration);
results.conditionTable=vertcat(pieces{:});
results.homoscedasticReferenceTable=homoscedasticTable;
if cfg.CreateFigures
    results.figureHandles=plotPairedProbePoissonSweep(results);
else
    results.figureHandles=gobjects(0);
end

if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_PoissonVariance_OffsetSummary.csv'));
    writetable(results.poolCalibrationTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_PoissonVariance_PoolCalibration.csv'));
    if ~isempty(homoscedasticTable)
        writetable(homoscedasticTable,fullfile(cfg.OutputDir, ...
            'PairedProbe_PoissonVariance_HomoscedasticReference.csv'));
    end
    returnedHandles=results.figureHandles; results.figureHandles=[];
    save(fullfile(cfg.OutputDir,'PairedProbe_PoissonVariance_Results.mat'), ...
        'results','-v7.3');
    results.figureHandles=returnedHandles;
    if cfg.CreateFigures
        savePairedProbePoissonFigures(results.figureHandles,cfg.OutputDir);
    end
end
end

function cfg=parseConfig(cfg,varargin)
if numel(varargin)==1 && isstruct(varargin{1})
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
validateattributes(cfg.ReferenceBaseDriftCoherence,{'numeric'},{'scalar','nonnegative'});
validateattributes(cfg.RandomMotionRateSpPerSec,{'numeric'},{'scalar','positive'});
validateattributes(cfg.AlignedRateSlopeSpPerSecPerPct,{'numeric'},{'scalar','positive'});
validateattributes(cfg.PreStepDurationSec,{'numeric'},{'scalar','positive'});
validateattributes(cfg.StepDurationSec,{'numeric'},{'scalar','positive'});
validateattributes(cfg.FanoFactor,{'numeric'},{'scalar','positive'});
end

function T=makeConditionTable(sim,threshold,cfg)
r=sim.rateDiagnostics; g=sim.gain; w=sim.winners; pcfg=cfg.PairedConfig;
T=table(sim.baseDriftCoherence,sim.offsetDeg,string(sim.yoking), ...
    sim.performance,sim.stepCoherence,threshold,sim.effectivePoolSize, ...
    g.drift,g.driftSE,g.probePerStream,g.probePerStreamSE, ...
    g.probeTotal,g.probeTotalSE,g.probePerStream/g.drift, ...
    w.changedDrift,w.changedProbe,w.changedOther,w.noChangeDrift, ...
    w.noChangeProbe,w.noChangeOther,sim.fit.intercept,sim.fit.logLikelihood, ...
    sim.fit.converged,r.preRateDrift,r.preRateProbe,r.preRateMinimum, ...
    r.preRateMaximum,r.preClampFraction,r.changedClampFraction, ...
    r.noChangeClampFraction,cfg.RandomMotionRateSpPerSec, ...
    cfg.AlignedRateSlopeSpPerSecPerPct,pcfg.InternalCorrelationLengthDeg, ...
    circularNoiseEffectiveRank(round(360/pcfg.BankSpacingDeg),pcfg), ...
    'VariableNames',{'BaseDriftCoherence','OffsetDeg','Yoking', ...
    'Performance','StepCoherence','Threshold75','EffectivePoolSize', ...
    'DriftGain','DriftGainSE','ProbeGainPerStream','ProbeGainPerStreamSE', ...
    'ProbeGainTotal','ProbeGainTotalSE','NormalizedProbePerStream', ...
    'ChangedDriftWinner','ChangedProbeWinner','ChangedOtherWinner', ...
    'NoChangeDriftWinner','NoChangeProbeWinner','NoChangeOtherWinner', ...
    'ProbitIntercept','LogLikelihood','ProbitConverged','PreRateDrift', ...
    'PreRateProbe','PreRateMinimum','PreRateMaximum','PreClampFraction', ...
    'ChangedClampFraction','NoChangeClampFraction','RandomMotionRate', ...
    'AlignedRateSlope','InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end
