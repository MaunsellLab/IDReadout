function results = runPairedProbeCorrelationSweep(varargin)
% RUNPAIREDPROBECORRELATIONSWEEP  Sweep smooth mechanism-noise correlation.
%
% COMPLETE MODEL SPECIFICATION FOR THE CORRELATION-LENGTH STAGE
% ============================================================
% Purpose
% -------
% The dense direction grid is a numerical sampling of a continuous response
% field, not a set of independent Gaussian competitors.  This stage measures
% how the smoothness and effective dimensionality of that field affect the
% paired-probe gain curves, winner identities, threshold invariance, and the
% internal-noise SD required to place a fixed 20% increment at 75% correct.
%
% Inherited delta-coherence model
% --------------------------------
% Every direction mechanism receives its tuning-weighted mean step coherence
% minus mean prestep coherence before hard-max pooling.  Constant drift and
% probe pedestals therefore cancel from the mechanism means.  The imposed
% drift and probe coherence changes retain their fixed scalar Gaussian noise.
% There is no downstream noise or lapse process.  See the complete header of
% runPairedProbeOffsetSweep.m for the inherited equations, gain recovery,
% winner classification, and fixed-winner/winner-switch decomposition.
%
% Quantity being swept
% --------------------
% Internal mechanism noise is a circular Gaussian field with nominal covariance
%
%   rho(d) = exp{kappa[cos(d)-1]},
%
% where kappa is set by InternalCorrelationLengthDeg.  The circular mean is
% removed on every trial and the marginal SD is restored.  This removal is held
% fixed throughout the present sweep; consequently, long-range correlations
% become negative.  A later control can examine mean removal separately.
%
% CorrelationLengthsDeg defaults to [10 15 20 30 45 60 90] deg.  Shorter
% lengths create more independent noise hills around the direction circle;
% longer lengths create fewer, broader hills.  The implemented covariance
% participation ratio, after circular-mean removal, is reported as
% EffectiveNoiseRank.  It is a dimensionality diagnostic, not an assertion
% that max statistics equal those of that many independent Gaussian variables.
%
% Calibration and common random numbers
% -------------------------------------
% At EACH correlation length, InternalNoiseSD is recalibrated once so the fixed
% drift increment (20% by default) produces 75% correct at the single reference
% offset/yoking condition.  That SD is then fixed across all offsets and both
% yoking modes for that length.  This separates changes in noise structure from
% trivial changes in the behavioral operating point.  Identical seeds generate
% common external stream draws and common underlying white Gaussian arrays
% across lengths; only their circular filtering changes.
%
% Default scope and outputs
% -------------------------
% To keep this exploratory sweep faster than the primary full-offset run, the
% defaults use 12,000 analysis trials, 8,000 calibration trials, and offsets
% [10 25 45 60 75 90 120 135 170 180].  Outputs combine every length into:
% an offset-summary table, decomposition table, noise-calibration table, MAT
% result, and figures for gain, operating point, calibrated SD/effective rank,
% selected offsets, and winner fractions.
%
% Usage
% -----
%   results = runPairedProbeCorrelationSweep;
%   results = runPairedProbeCorrelationSweep( ...
%       'CorrelationLengthsDeg',[15 30 60]);
%   cfg = pairedProbeCorrelationSweepDefaults;
%   cfg.PairedConfig.NTrials = 20000;
%   results = runPairedProbeCorrelationSweep(cfg);

cfg=pairedProbeCorrelationSweepDefaults();
cfg=parseSweepConfig(cfg,varargin{:});
if cfg.SaveOutputs && ~exist(cfg.OutputDir,'dir'), mkdir(cfg.OutputDir); end

nLengths=numel(cfg.CorrelationLengthsDeg);
conditionPieces=cell(nLengths,1);
decompositionPieces=cell(nLengths,1);
calibrationPieces=cell(nLengths,1);

for ll=1:nLengths
    stageCfg=cfg.PairedConfig;
    stageCfg.InternalCorrelationLengthDeg=cfg.CorrelationLengthsDeg(ll);
    stageCfg.Visible='off';
    stageCfg.CreateFigures=false;
    stageCfg.SaveOutputs=false;
    stage=runPairedProbeOffsetSweep(stageCfg);
    conditionPieces{ll}=stage.conditionTable;
    decompositionPieces{ll}=stage.decompositionTable;
    calibrationPieces{ll}=stage.noiseCalibrationTable;
end

results=struct();
results.config=cfg;
results.conditionTable=vertcat(conditionPieces{:});
results.decompositionTable=vertcat(decompositionPieces{:});
results.noiseCalibrationTable=vertcat(calibrationPieces{:});
results.figureHandles=plotPairedProbeCorrelationSweep(results);

if cfg.SaveOutputs
    writetable(results.conditionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_CorrelationSweep_OffsetSummary.csv'));
    writetable(results.decompositionTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_CorrelationSweep_Decomposition.csv'));
    writetable(results.noiseCalibrationTable,fullfile(cfg.OutputDir, ...
        'PairedProbe_CorrelationSweep_NoiseCalibration.csv'));
    returnedFigureHandles=results.figureHandles;
    results.figureHandles=[];
    save(fullfile(cfg.OutputDir,'PairedProbe_CorrelationSweep_Results.mat'), ...
        'results','-v7.3');
    results.figureHandles=returnedFigureHandles;
    savePairedProbeCorrelationSweepFigures(results.figureHandles,cfg.OutputDir);
end
end

function cfg=parseSweepConfig(cfg,varargin)
if isempty(varargin)
    % Use defaults.
elseif numel(varargin)==1 && isstruct(varargin{1})
    supplied=varargin{1}; names=fieldnames(supplied);
    for k=1:numel(names)
        if ~isfield(cfg,names{k}), error('Unknown sweep field: %s',names{k}); end
        cfg.(names{k})=supplied.(names{k});
    end
elseif mod(numel(varargin),2)==0
    for k=1:2:numel(varargin)
        name=char(varargin{k});
        if ~isfield(cfg,name), error('Unknown sweep field: %s',name); end
        cfg.(name)=varargin{k+1};
    end
else
    error('Use a sweep configuration struct or name/value pairs.');
end
validateattributes(cfg.CorrelationLengthsDeg,{'numeric'}, ...
    {'vector','positive','finite'});
if ~isstruct(cfg.PairedConfig)
    error('PairedConfig must be a structure returned by pairedProbeDefaults.');
end
end
