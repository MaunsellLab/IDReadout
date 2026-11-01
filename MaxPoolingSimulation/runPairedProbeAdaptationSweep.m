function results = runPairedProbeAdaptationSweep(varargin)
% RUNPAIREDPROBEADAPTATIONSWEEP  Drift cancellation plus selective adaptation.
%
% COMPLETE MODEL SPECIFICATION FOR THE ADAPTATION STAGE
% =====================================================
% Scientific question
% -------------------
% Exact mechanism-wise delta coherence stabilizes threshold across base-drift
% coherence, but removes the probe pedestals that generate directional
% convergence.  Absolute-response temporal max preserves those pedestals but
% lets a 40% base drift monopolize the max.  This stage tests a biologically
% constrained alternative: cancel the perceived sustained drift exactly, keep
% the probe pedestals, and add firing-rate-dependent adaptation centered on the
% prestep drift.  It replaces, rather than supplements, the absolute-response
% temporal-max model.
%
% Mechanism response
% ------------------
% Preferred directions phi_m densely tile the circle (1-deg default).  Drift
% is at 0 deg and paired probes are at +/-offset.  Each stream reaches every
% mechanism through centered wrapped-Gaussian MT sensitivity H with sigma
% 37.5 deg.  In either spatial patch the step-period landscape is
%
%   E_m = sum_{j=probes} pedestal_j H(theta_j-phi_m)
%         + sum_{j=drift,+probe,-probe} deltaC_j H(theta_j-phi_m)
%         + step*H(-phi_m) - A(base) G(phi_m) + epsilon_m .
%
% The fixed positive drift step is present only in the changed patch.  Exact
% perceptual cancellation removes (BaseDriftCoherence+DriftPedestal)H from
% both patches before the max; it does not remove the two probe pedestals.
% Thus sustained drift cannot privilege the drift mechanism merely by being
% large, while probe convergence remains available.
%
% Adaptation is a nonnegative suppressive field
%
%   G(phi)=exp[-distance(phi,0)^2/(2*sigma_A^2)]
%   A(base)=AftereffectAt60*(BaseDriftCoherence+DriftPedestal)/60 .
%
% sigma_A equals the MT tuning sigma by default.  AftereffectAt60 is expressed
% in coherence-equivalent percentage points at a 60% prestep.  The default
% sweep 0:8 samples the moderate range implicated by the first adaptation run
% and brackets much of the approximately 10%-coherence horizontal shift
% suggested by the preferred/null MT lines after a 60% prestep.  Mechanisms
% far from drift are spared, not actively excited; any relative enhancement
% occurs naturally through max competition.  Linear base scaling is the first
% firing-rate-dependent approximation and is explicit for later replacement.
%
% Noise, max, and choice
% ----------------------
% One scalar Gaussian epoch-mean fluctuation is drawn per drift/probe stream
% and patch.  Probe streams are independent or yoked; the two spatial patches
% remain independent.  Additive mechanism noise is a smooth circular Gaussian
% field with 30-deg nominal correlation length and circular-mean removal.  It
% is not independent noise at every grid direction.  There is no downstream
% noise and no lapse process.  The abstraction has no frame-resolved coherence
% and does not enforce a fixed dot budget.
%
% A hard direction max is computed separately in changed and no-change
% patches; the model chooses changed when its maximum is larger.  Soft max and
% an eligibility field are intentionally deferred so this stage isolates the
% effect of measured adaptation.
%
% Calibration and predictions
% ---------------------------
% One internal-noise SD is calibrated so a fixed 15.3% increment gives 75%
% correct at base 7.071%, offset 120 deg, independent probes, and zero
% adaptation.  The 15.3% anchor is the current IDS mean from 17 sessions and
% is substantially more reliable than the 11.2% 170-deg mean from 3 sessions.
% The calibrated SD is held fixed across both base drifts,
% every offset, yoking mode, and adaptation strength.  Every condition reports
% performance at the common 15.3% step plus a separately measured 75% threshold.
% Threshold stability and gain structure are therefore predictions outside
% the single calibration condition.  The output also reports the matched
% 40%-minus-7.071% threshold difference.  Approximately 5 coherence points is
% the provisional empirical target based on IDS versus IDR threshold means;
% it is not assumed to be a pure base-coherence effect because the experiments
% also differ in yoking, direction predictability, and other task details.
%
% Gain and winner diagnostics
% ---------------------------
% Correct/error choices are fit by simultaneous probit regression on the
% changed-minus-no-change scalar stream fluctuations.  For independent probes,
% the two coefficients are averaged per stream; for yoked probes, half of the
% identifiable pair coefficient is reported.  Counterfactual perturbations
% divide derivatives into fixed-winner and winner-switch terms.  Winners are
% also reported as drift/probe/other categories and in angular bins with local
% drift and probe sensitivities and signed probability-weighted contributions.
% The angular fixed-winner totals omit switch terms and decision-boundary
% weighting, so they need not equal recovered probit gains.
%
% Default sweep and usage
% -----------------------
% Bases: 7.071 and 40%; offsets: the five current IDS values 10, 20, 45,
% 120, and 170 deg; independent and yoked probes; aftereffects at 60%: 0:8%.
% The operating-point figures overlay IDS threshold means [13.0 11.9 16.2
% 15.3 11.2] and the provisional 21-22% IDR range.  For example:
%
%   results = runPairedProbeAdaptationSweep;
%   cfg = adaptationSweepDefaults;
%   cfg.AftereffectAt60Values = [0 5 10];
%   results = runPairedProbeAdaptationSweep(cfg);

cfg=adaptationSweepDefaults(); cfg=parseConfig(cfg,varargin{:});
pcfg=cfg.PairedConfig;
pcfg.EffectiveNoiseRank=circularNoiseEffectiveRank(round(360/pcfg.BankSpacingDeg),pcfg);
cfg.PairedConfig=pcfg;
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

calibrationBank=makeDensePairedTrialBank(pcfg.NCalibrationTrials,pcfg,pcfg.Seed+610000);
analysisBank=makeDensePairedTrialBank(pcfg.NTrials,pcfg,pcfg.Seed+620000);
if pcfg.CalibrateInternalNoise
    calibration=calibrateAdaptationNoise(calibrationBank,cfg);
    pcfg.InternalNoiseSD=calibration.internalNoiseSD; cfg.PairedConfig=pcfg;
else
    calibration=struct('internalNoiseSD',pcfg.InternalNoiseSD,'performance',NaN, ...
        'performanceWithNoInternalNoise',NaN,'baseDriftCoherence', ...
        cfg.ReferenceBaseDriftCoherence,'aftereffectAt60', ...
        cfg.ReferenceAftereffectAt60,'actualAftereffect',NaN, ...
        'offsetDeg',pcfg.NoiseCalibrationOffsetDeg, ...
        'yoked',logical(pcfg.NoiseCalibrationYoked), ...
        'stepCoherence',pcfg.FixedStepCoherence);
end
calibration.correlationLengthDeg=pcfg.InternalCorrelationLengthDeg;
calibration.effectiveNoiseRank=pcfg.EffectiveNoiseRank;
calibration.removeCircularMean=logical(pcfg.RemoveInternalCircularMean);

n=numel(cfg.BaseDriftCoherences)*numel(pcfg.OffsetsDeg)*2*numel(cfg.AftereffectAt60Values);
C=cell(n,1); D=cell(2*n,1); W=cell(n,1); piece=0; dd=0;
for base=cfg.BaseDriftCoherences
    for offset=pcfg.OffsetsDeg
        for yoked=[false true]
            for ae=cfg.AftereffectAt60Values
                piece=piece+1;
                threshold=calibrateAdaptationStep(calibrationBank,offset,yoked,base,ae, ...
                    pcfg.InternalNoiseSD,cfg);
                sim=simulateAdaptationCondition(analysisBank,offset,yoked, ...
                    pcfg.FixedStepCoherence,base,ae,pcfg.InternalNoiseSD,cfg);
                C{piece}=conditionTable(sim,threshold,cfg);
                W{piece}=winnerDirectionTable(sim,cfg);
                for stream=1:2
                    dd=dd+1; D{dd}=decompositionTable(sim,stream,cfg);
                end
            end
        end
    end
end
results=struct('config',cfg,'noiseCalibration',calibration, ...
    'noiseCalibrationTable',struct2table(calibration), ...
    'conditionTable',vertcat(C{:}),'decompositionTable',vertcat(D{:}), ...
    'winnerDirectionTable',vertcat(W{:}));
results.baseContrastTable=buildBaseContrastTable(results.conditionTable);
if cfg.CreateFigures, results.figureHandles=plotPairedProbeAdaptationSweep(results);
else, results.figureHandles=gobjects(0); end
if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir,'PairedProbe_Adaptation_ConditionSummary.csv'));
    writetable(results.decompositionTable,fullfile(cfg.OutputDir,'PairedProbe_Adaptation_Decomposition.csv'));
    writetable(results.winnerDirectionTable,fullfile(cfg.OutputDir,'PairedProbe_Adaptation_WinnerDirectionContributions.csv'));
    writetable(results.noiseCalibrationTable,fullfile(cfg.OutputDir,'PairedProbe_Adaptation_NoiseCalibration.csv'));
    writetable(results.baseContrastTable,fullfile(cfg.OutputDir,'PairedProbe_Adaptation_BaseContrast.csv'));
    handles=results.figureHandles; results.figureHandles=[];
    save(fullfile(cfg.OutputDir,'PairedProbe_Adaptation_Results.mat'),'results','-v7.3');
    results.figureHandles=handles;
    if cfg.CreateFigures, savePairedProbeAdaptationFigures(handles,cfg.OutputDir); end
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
else, error('Use a configuration struct or name/value pairs.');
end
pcfg=cfg.PairedConfig;
validateattributes(cfg.BaseDriftCoherences,{'numeric'},{'vector','nonnegative'});
validateattributes(cfg.AftereffectAt60Values,{'numeric'},{'vector','nonnegative'});
validateattributes(cfg.AdaptationReferenceCoherence,{'numeric'},{'scalar','positive'});
validateattributes(cfg.AdaptationWidthDeg,{'numeric'},{'scalar','positive','<=',180});
validateattributes(cfg.EmpiricalIDSOffsetsDeg,{'numeric'},{'vector','>=',0,'<=',180});
validateattributes(cfg.EmpiricalIDSThresholds,{'numeric'},{'vector','positive'});
if numel(cfg.EmpiricalIDSOffsetsDeg)~=numel(cfg.EmpiricalIDSThresholds)
    error('Empirical IDS offsets and thresholds must have equal lengths.');
end
if numel(cfg.EmpiricalIDSSessionCounts)~=numel(cfg.EmpiricalIDSThresholds)
    error('Empirical IDS session counts must match the threshold vector.');
end
validateattributes(cfg.EmpiricalIDRThresholdRange,{'numeric'},{'numel',2,'positive'});
validateattributes(cfg.ExpectedBaseThresholdIncrease,{'numeric'},{'scalar','nonnegative'});
validateattributes(cfg.WinnerDirectionBinWidthDeg,{'numeric'},{'scalar','positive','<=',360});
if abs(360/cfg.WinnerDirectionBinWidthDeg-round(360/cfg.WinnerDirectionBinWidthDeg))>1e-10
    error('WinnerDirectionBinWidthDeg must divide 360 exactly.');
end

validateattributes(pcfg.OffsetsDeg,{'numeric'},{'vector','>=',0,'<=',180});
if abs(360/pcfg.BankSpacingDeg-round(360/pcfg.BankSpacingDeg))>1e-10
    error('BankSpacingDeg must divide 360 exactly.');
end
end

function B=buildBaseContrastTable(T)
bases=sort(unique(T.BaseDriftCoherence));
if numel(bases)~=2
    error('BaseContrastTable currently requires exactly two base coherences.');
end
low=T(T.BaseDriftCoherence==bases(1),:);
high=T(T.BaseDriftCoherence==bases(2),:);
keys={'AftereffectAt60','OffsetDeg','Yoking'};
low=sortrows(low,keys); high=sortrows(high,keys);
if height(low)~=height(high) || any(low.AftereffectAt60~=high.AftereffectAt60) || ...
        any(low.OffsetDeg~=high.OffsetDeg) || any(low.Yoking~=high.Yoking)
    error('Low- and high-base conditions do not align.');
end
B=table(repmat(bases(1),height(low),1),repmat(bases(2),height(low),1), ...
    low.AftereffectAt60,low.OffsetDeg,low.Yoking,low.Threshold75, ...
    high.Threshold75,high.Threshold75-low.Threshold75,low.Performance, ...
    high.Performance,low.DriftGain,high.DriftGain,low.ProbeGainPerStream, ...
    high.ProbeGainPerStream,low.NormalizedProbePerStream, ...
    high.NormalizedProbePerStream,'VariableNames', ...
    {'LowBaseDriftCoherence','HighBaseDriftCoherence','AftereffectAt60', ...
    'OffsetDeg','Yoking','LowBaseThreshold75','HighBaseThreshold75', ...
    'ThresholdIncrease','LowBasePerformance','HighBasePerformance', ...
    'LowBaseDriftGain','HighBaseDriftGain','LowBaseProbeGainPerStream', ...
    'HighBaseProbeGainPerStream','LowBaseNormalizedProbePerStream', ...
    'HighBaseNormalizedProbePerStream'});
end

function T=conditionTable(s,threshold,cfg)
g=s.gain; w=s.winners; p=cfg.PairedConfig;
T=table(s.baseDriftCoherence,s.aftereffectAt60,s.actualAftereffect, ...
    s.cancelledDriftCoherence,s.offsetDeg,string(s.yoking),s.performance, ...
    s.stepCoherence,threshold,g.drift,g.driftSE,g.probePerStream, ...
    g.probePerStreamSE,g.probeTotal,g.probeTotalSE,g.probePerStream/g.drift, ...
    w.changedDrift,w.changedProbe,w.changedOther,w.noChangeDrift, ...
    w.noChangeProbe,w.noChangeOther,s.fit.intercept,s.fit.logLikelihood, ...
    s.fit.converged,p.InternalNoiseSD,p.ProbePedestal, ...
    cfg.AdaptationReferenceCoherence,cfg.AdaptationWidthDeg, ...
    p.InternalCorrelationLengthDeg,p.EffectiveNoiseRank, ...
    'VariableNames',{'BaseDriftCoherence','AftereffectAt60','ActualAftereffect', ...
    'CancelledDriftCoherence','OffsetDeg','Yoking','Performance','StepCoherence', ...
    'Threshold75','DriftGain','DriftGainSE','ProbeGainPerStream', ...
    'ProbeGainPerStreamSE','ProbeGainTotal','ProbeGainTotalSE', ...
    'NormalizedProbePerStream','ChangedDriftWinner','ChangedProbeWinner', ...
    'ChangedOtherWinner','NoChangeDriftWinner','NoChangeProbeWinner', ...
    'NoChangeOtherWinner','ProbitIntercept','LogLikelihood','ProbitConverged', ...
    'InternalNoiseSD','ProbePedestal','AdaptationReferenceCoherence', ...
    'AdaptationWidthDeg','InternalCorrelationLengthDeg','EffectiveNoiseRank'});
end

function T=decompositionTable(s,stream,cfg)
if stream==1, name="Drift"; d=s.decomposition.drift;
else, name="ProbePerStream"; d=s.decomposition.probePerStream; end
T=table(s.baseDriftCoherence,s.aftereffectAt60,s.actualAftereffect,s.offsetDeg, ...
    string(s.yoking),name,d.changedFixed,d.changedSwitch,d.noChangeFixed, ...
    d.noChangeSwitch,d.decisionFixed,d.decisionSwitch,d.decisionTotal, ...
    d.changedSwitchProbability,d.noChangeSwitchProbability, ...
    cfg.PairedConfig.InternalNoiseSD,'VariableNames', ...
    {'BaseDriftCoherence','AftereffectAt60','ActualAftereffect','OffsetDeg', ...
    'Yoking','Stream','ChangedFixedWinner','ChangedWinnerSwitch', ...
    'NoChangeFixedWinner','NoChangeWinnerSwitch','DecisionFixedWinner', ...
    'DecisionWinnerSwitch','DecisionTotal','ChangedSwitchProbability', ...
    'NoChangeSwitchProbability','InternalNoiseSD'});
end

function T=winnerDirectionTable(s,cfg)
d=s.winnerDirection; n=numel(d.directionDeg);
T=table(repmat(s.baseDriftCoherence,n,1),repmat(s.aftereffectAt60,n,1), ...
    repmat(s.actualAftereffect,n,1),repmat(s.offsetDeg,n,1), ...
    repmat(string(s.yoking),n,1),d.directionDeg,d.changedWinnerFraction, ...
    d.noChangeWinnerFraction,d.patchAverageWinnerFraction, ...
    d.driftLocalSensitivity,d.probeLocalSensitivityPerStream, ...
    d.driftWinnerWeightedContribution,d.probeWinnerWeightedContributionPerStream, ...
    repmat(d.driftFixedWinnerTotal,n,1), ...
    repmat(d.probeFixedWinnerTotalPerStream,n,1), ...
    repmat(cfg.WinnerDirectionBinWidthDeg,n,1),'VariableNames', ...
    {'BaseDriftCoherence','AftereffectAt60','ActualAftereffect','OffsetDeg', ...
    'Yoking','WinnerDirectionDeg','ChangedWinnerFraction', ...
    'NoChangeWinnerFraction','PatchAverageWinnerFraction','DriftLocalSensitivity', ...
    'ProbeLocalSensitivityPerStream','DriftWinnerWeightedContribution', ...
    'ProbeWinnerWeightedContributionPerStream','DriftFixedWinnerTotal', ...
    'ProbeFixedWinnerTotalPerStream','DirectionBinWidthDeg'});
end
