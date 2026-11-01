function results = runPairedProbeOffsetSweep(varargin)
% RUNPAIREDPROBEOFFSETSWEEP  Change-detection max pooling with paired probes.
%
% COMPLETE MODEL SPECIFICATION FOR THE DELTA-COHERENCE PAIRED-PROBE STAGE
% =======================================================================
% Scientific purpose
% ------------------
% This stage tests a change-detection max-pooling model.  It separates a dense
% bank of eligible direction mechanisms from three experimentally imposed
% stimulus streams: drift at 0 deg and paired probes at +/-offset.  It asks how
% tuning overlap, probe-noise yoking, and max competition affect recovered
% drift and probe gains after every mechanism has removed its own prestep
% baseline.
%
% Prestep-to-step coherence change
% --------------------------------
% Each mechanism is fed the difference between its tuning-weighted mean
% coherence in the step and prestep epochs.  Differencing occurs independently
% within every direction mechanism BEFORE max pooling:
%
%   Delta C_m = sum_j H(theta_j-phi_m) * (C_step,j-C_pre,j).
%
% Temporal frames and detailed prestep fluctuations are not simulated.  The
% three stream directions are theta = [0,+delta,-delta].  Their prestep means
% are BaseDriftCoherence+DriftPedestal for the drift stream and ProbePedestal
% for each probe stream.  Step means contain the same constants plus one scalar
% Gaussian coherence change x_j per stream.  The changed patch additionally
% receives a fixed positive drift increment S; the no-change patch does not.
% Consequently, the constant drift base and all pedestals cancel exactly.
% BaseDriftCoherence, DriftPedestal, and ProbePedestal remain configuration
% fields so this invariance can be verified explicitly, but they cannot
% privilege a max competitor in this linear delta-coherence stage.
%
% S defaults to 20 coherence-percentage units.  Drift change-noise SD defaults
% to 10/sqrt(5); probe change-noise SD defaults to 7.071/sqrt(5).  These scalar
% Gaussian variables represent only the effective difference between epoch
% means.  They do not model 50-ms frames, temporal ordering, or separate noisy
% estimates of the two epoch means.  The abstraction does not enforce a 100%
% dot budget.
%
% In the Independent condition x_+ and x_- are independent.  In the Yoked
% condition x_+ = x_- within each patch.  Changed and no-change patches remain
% independent.  At 180 deg the two probe directions coincide; their coherence
% changes and response contributions add.
%
% Dense MT response bank
% ----------------------
% Eligible preferred directions phi span the full circle at
% Config.BankSpacingDeg (1 deg by default).  Each stream affects every
% mechanism through
%
%   H(theta-phi) = exp[-(theta-phi)^2/(2*sigma_T^2)] - mean_circle(T),
%
% with wrapped angular distance and sigma_T = 37.5 deg by default.  There is
% no readout-weight or eligibility difference across directions.  All
% directions remain eligible and can win through tuning overlap, intermediate
% response peaks, stimulus-change noise, or internal noise.
%
% Internal mechanism noise
% ------------------------
% Additive Gaussian internal noise is a smooth circular field, not 360
% independent samples.  Its covariance is a positive-definite von-Mises form
%
%   Corr(epsilon(phi),epsilon(phi+d)) = exp{kappa[cos(d)-1]},
%
% where kappa is chosen from Config.InternalCorrelationLengthDeg.  The circular
% mean is removed on every trial to prevent a patch-common offset from acting
% as downstream decision noise.  This removal also creates negative long-range
% correlations; EffectiveNoiseRank reports the participation ratio of the
% covariance after that operation.  There is no lapse process or downstream
% decision noise.
%
% Fixed operating point and one-time noise calibration
% ----------------------------------------------------
% Each patch value is the hard maximum over the dense response bank.  Choice
% is correct when the changed-patch maximum exceeds the no-change maximum.
% The internal-noise SD is calibrated ONCE: FixedStepCoherence (20% by default)
% is matched to TargetPerformance (75%) at NoiseCalibrationOffsetDeg and
% NoiseCalibrationYoked.  External stream-change noise remains fixed.  The
% resulting internal-noise SD is then held fixed across every offset and both
% yoking conditions.  If CalibrateInternalNoise is false, the supplied
% InternalNoiseSD is used without calibration.
%
% Every condition is analyzed at the same FixedStepCoherence, so reported
% performance is a prediction rather than a forced operating point.  A separate
% Threshold75 is measured for each condition with the already-fixed noise; this
% diagnoses whether the model produces the observed stable approximately 20%
% increment threshold.
%
% Correct/error choices are fit with a simultaneous probit regression on
% change-minus-no-change stream deltas.  Independent probes receive separate
% coefficients, summarized by their mean and sum.  The yoked common predictor
% receives one identifiable coefficient; this is the total pair effect, and
% half is reported as the per-stream equivalent.  Gains are recovered from
% choices rather than read from tuning profiles.
%
% Fixed-winner versus winner-switch decomposition
% ------------------------------------------------
% A small paired counterfactual perturbation is applied to each stream profile
% while all trial noise is held fixed.  For each patch the change in its maximum
% is divided into a fixed-winner effect and a winner-switch effect.  Changed-
% and no-change-patch contributions are retained separately.  For yoked probes,
% the perturbation changes both probe streams together and half of the pair
% value is reported as the per-stream equivalent.
%
% Winner-location diagnostic
% --------------------------
% Winning directions are classified as drift-associated, probe-associated, or
% other/intermediate using Config.WinnerAssociationRadiusDeg.  This shows
% whether probe gain is carried by probe-centered winners or indirectly through
% intermediate mechanisms.  Output tables retain the historical label Other.
%
% Later extensions
% ----------------
% Later stages may replace the scalar epoch-mean changes with empirical
% leave-one-out kernel-weighted predictors, allow imperfect baseline removal,
% or introduce response-dependent noise.  None is included here.
%
% Usage
% -----
%   results = runPairedProbeOffsetSweep;
%   results = runPairedProbeOffsetSweep('FixedStepCoherence',20, ...
%       'NoiseCalibrationOffsetDeg',90,'Visible','on');
%   cfg = pairedProbeDefaults; cfg.OffsetsDeg = [10 25 45 120 170 179];
%   results = runPairedProbeOffsetSweep(cfg);

cfg = pairedProbeDefaults();
cfg = parsePairedConfig(cfg,varargin{:});
cfg.EffectiveNoiseRank=circularNoiseEffectiveRank( ...
    round(360/cfg.BankSpacingDeg),cfg);
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

calibrationBank = makeDensePairedTrialBank(cfg.NCalibrationTrials,cfg, ...
    cfg.Seed+100000);
analysisBank = makeDensePairedTrialBank(cfg.NTrials,cfg,cfg.Seed+200000);

if cfg.CalibrateInternalNoise
    noiseCalibration=calibratePairedProbeNoise(calibrationBank,cfg);
    cfg.InternalNoiseSD=noiseCalibration.internalNoiseSD;
else
    noiseCalibration=struct('internalNoiseSD',cfg.InternalNoiseSD, ...
        'performance',NaN,'performanceWithNoInternalNoise',NaN, ...
        'offsetDeg',cfg.NoiseCalibrationOffsetDeg, ...
        'yoked',logical(cfg.NoiseCalibrationYoked), ...
        'stepCoherence',cfg.FixedStepCoherence);
end
noiseCalibration.correlationLengthDeg=cfg.InternalCorrelationLengthDeg;
noiseCalibration.effectiveNoiseRank=cfg.EffectiveNoiseRank;
noiseCalibration.removeCircularMean=logical(cfg.RemoveInternalCircularMean);

nConditions = numel(cfg.OffsetsDeg)*2;
conditionPieces = cell(nConditions,1);
decompositionPieces = cell(2*nConditions,1);
piece = 0; decompPiece = 0;

for oo = 1:numel(cfg.OffsetsDeg)
    offset = cfg.OffsetsDeg(oo);
    for yoked = [false true]
        piece = piece+1;
        threshold75 = calibratePairedProbeStep(calibrationBank,offset,yoked,cfg);
        sim = simulatePairedProbeCondition(analysisBank,offset,yoked, ...
            cfg.FixedStepCoherence,cfg);
        sim.threshold75=threshold75;
        conditionPieces{piece} = conditionTable(sim,cfg);
        for ss = 1:2
            decompPiece = decompPiece+1;
            decompositionPieces{decompPiece} = decompositionTable(sim,ss,cfg);
        end
    end
end

results = struct();
results.config = cfg;
results.noiseCalibration = noiseCalibration;
results.noiseCalibrationTable = struct2table(noiseCalibration);
results.conditionTable = vertcat(conditionPieces{:});
results.decompositionTable = vertcat(decompositionPieces{:});
if cfg.CreateFigures
    results.figureHandles = plotPairedProbeSweep(results);
else
    results.figureHandles = gobjects(0);
end

if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_OffsetSummary.csv'));
    writetable(results.decompositionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_Decomposition.csv'));
    writetable(results.noiseCalibrationTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_NoiseCalibration.csv'));
    returnedFigureHandles = results.figureHandles;
    results.figureHandles = [];
    save(fullfile(cfg.OutputDir,'PairedProbe_Results.mat'),'results','-v7.3');
    results.figureHandles = returnedFigureHandles;
    if cfg.CreateFigures
        savePairedProbeFigures(results.figureHandles,cfg.OutputDir);
    end
end
end

function cfg = parsePairedConfig(cfg,varargin)
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
validateattributes(cfg.OffsetsDeg,{'numeric'},{'vector','>=',0,'<=',180});
validateattributes(cfg.BankSpacingDeg,{'numeric'},{'scalar','positive','<=',10});
validateattributes(cfg.FixedStepCoherence,{'numeric'},{'scalar','positive'});
validateattributes(cfg.NoiseCalibrationOffsetDeg,{'numeric'}, ...
    {'scalar','>=',0,'<=',180});
validateattributes(cfg.InternalNoiseSD,{'numeric'},{'scalar','nonnegative'});
validateattributes(cfg.NoiseCalibrationYoked,{'logical','numeric'},{'scalar'});
if abs(360/cfg.BankSpacingDeg-round(360/cfg.BankSpacingDeg))>1e-10
    error('BankSpacingDeg must divide 360 exactly.');
end
end

function T = conditionTable(sim,cfg)
T = table(sim.offsetDeg,string(sim.yoking),sim.performance,sim.stepCoherence, ...
    sim.threshold75,sim.gain.drift,sim.gain.driftSE,sim.gain.probePerStream, ...
    sim.gain.probePerStreamSE,sim.gain.probeTotal,sim.gain.probeTotalSE, ...
    sim.gain.probePerStream/sim.gain.drift,sim.gain.probeTotal/sim.gain.drift, ...
    sim.winners.changedDrift,sim.winners.changedProbe,sim.winners.changedOther, ...
    sim.winners.noChangeDrift,sim.winners.noChangeProbe,sim.winners.noChangeOther, ...
    sim.fit.intercept,sim.fit.logLikelihood,sim.fit.converged, ...
    cfg.InternalNoiseSD,cfg.BaseDriftCoherence,cfg.DriftPedestal,cfg.ProbePedestal, ...
    cfg.InternalCorrelationLengthDeg,cfg.EffectiveNoiseRank, ...
    'VariableNames',{'OffsetDeg','Yoking','Performance','StepCoherence', ...
    'Threshold75','DriftGain','DriftGainSE','ProbeGainPerStream','ProbeGainPerStreamSE', ...
    'ProbeGainTotal','ProbeGainTotalSE','NormalizedProbePerStream', ...
    'NormalizedProbeTotal','ChangedDriftWinner','ChangedProbeWinner', ...
    'ChangedOtherWinner','NoChangeDriftWinner','NoChangeProbeWinner', ...
    'NoChangeOtherWinner','ProbitIntercept','LogLikelihood','ProbitConverged', ...
    'InternalNoiseSD','BaseDriftCoherence','DriftPedestal','ProbePedestal', ...
    'InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end

function T = decompositionTable(sim,streamNumber,cfg)
if streamNumber==1
    name="Drift"; d=sim.decomposition.drift;
else
    name="ProbePerStream"; d=sim.decomposition.probePerStream;
end
T=table(sim.offsetDeg,string(sim.yoking),name,d.changedFixed,d.changedSwitch, ...
    d.noChangeFixed,d.noChangeSwitch,d.decisionFixed,d.decisionSwitch, ...
    d.decisionTotal,d.changedSwitchProbability,d.noChangeSwitchProbability, ...
    cfg.InternalCorrelationLengthDeg,cfg.EffectiveNoiseRank, ...
    'VariableNames',{'OffsetDeg','Yoking','Stream','ChangedFixedWinner', ...
    'ChangedWinnerSwitch','NoChangeFixedWinner','NoChangeWinnerSwitch', ...
    'DecisionFixedWinner','DecisionWinnerSwitch','DecisionTotal', ...
    'ChangedSwitchProbability','NoChangeSwitchProbability', ...
    'InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end
