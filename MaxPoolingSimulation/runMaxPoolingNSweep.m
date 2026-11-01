function results = runMaxPoolingNSweep(varargin)
% RUNMAXPOOLINGNSWEEP  Sparse direction/max-pooling education simulation.
%
% MODEL AND STAGED SIMULATION SPECIFICATION
% =========================================
% Purpose
% -------
% This simulation asks how regression-derived stimulus gains emerge from a
% max readout even when every stimulus stream has the same direct alignment
% with one eligible direction mechanism.  It is deliberately simpler than
% the later IDR/IDQ/IDS trial-level model, but it preserves the ingredients
% needed to study mechanism number, MT tuning overlap, internal noise,
% operating point, and proximity to the within-patch maximum.
%
% Stage 1: sparse, perfectly stream-aligned direction mechanisms (implemented)
% ---------------------------------------------------------------------------
% For each n in Config.NValues, n preferred directions are placed evenly
% around the circle:
%
%       phi_j = 360*(j-1)/n,   j = 1,...,n.
%
% The first mechanism is aligned with the drift direction (0 deg).  There is
% one zero-mean scalar stimulus-noise integral x_j centered exactly on every
% mechanism.  There is therefore NO probe/mechanism misalignment and NO
% averaging over rotations.  Increasing n intentionally decreases the
% angular spacing and increases overlap/encroachment among streams.
%
% Each stream affects every mechanism through a centered MT sensitivity:
%
%       H(delta) = exp(-delta^2/(2*sigma_T^2)) - mean_circle(T).
%
% For patch q and mechanism m, after integration over the coherence step,
%
%       R_qm = sum_j x_qj H(phi_j-phi_m) + epsilon_qm.
%
% epsilon_qm is independent Gaussian mechanism-response noise.  The changed
% patch additionally receives deterministic coherence S in the drift stream:
%
%       R_change,m <- R_change,m + S H(0-phi_m).
%
% The patch value is the hard maximum over mechanisms.  Choice is correct
% when max_m R_change,m > max_m R_nochange,m.  There is no downstream decision
% noise and no lapse process.  Config.StreamPedestal can add a scalar pedestal
% to every stream (or a vector with one value per stream).  A common scalar
% pedestal is intentionally an invariance control: because the aligned tuning
% matrix is circulant, it adds the same response constant to every mechanism
% and cannot affect choices or gains.
%
% Noise and operating point
% -------------------------
% Stimulus streams are independent Gaussian scalar integrals with SD
% Config.ExternalNoiseSD.  This SD defines the stimulus unit; it can later be
% replaced by the empirically expected kernel-integral SD.  Mechanism noise
% has SD Config.InternalNoiseSD.  Their ratio is explicit and can be swept.
%
% Two operating modes are simulated from common model rules:
%   FIXED   - one S is calibrated to TargetPerformance at ReferenceN and then
%             held fixed across n.  Performance is allowed to change with n.
%   MATCHED - S is separately calibrated at every n to TargetPerformance.
%             Internal and external noise scales remain fixed.
%
% Calibration uses independent trials from the trials used to recover gains.
% Common random numbers are used within each calibration so the hard-max
% performance function is stable during bisection.
%
% Gain recovery
% -------------
% Correct/error choices are regressed simultaneously on the change-minus-
% no-change scalar noise integral for every stream, using a probit GLM fitted
% by toolbox-independent IRLS.  The fitted stream coefficients are the
% measured gains.  Symmetric streams at equal circular offset are averaged;
% normalized gain is the averaged coefficient divided by drift-stream gain.
% These coefficients are outputs of simulated choices, not model parameters.
%
% Stage 2 diagnostic: proximity to the maximum (implemented)
% -----------------------------------------------------------
% For mechanism m on each changed-patch trial, the competitive margin is
%
%       gap_m = R_m - max_{k~=m} R_k.
%
% A paired counterfactual perturbation adds a small delta directly to R_m
% while holding every random draw fixed.  The resulting change in the patch
% maximum, divided by delta, is plotted against gap_m.  The simulation also
% relates each aligned mechanism's win probability to its recovered stream
% gain.  These diagnostics expose the close-to-max gate without conditioning
% a choice regression on the same noise that selected the winner.
%
% Later stages
% ------------
%   A. IMPLEMENTED SEPARATELY in runPairedProbeOffsetSweep: symmetric probes
%      at +/-offset, probe pedestals, independent versus yoked probe noise, a
%      dense direction field, smooth correlated internal noise, and causal
%      fixed-winner/winner-switch decompositions.
%   B. Still later: systematic control and recovery of the dense field's
%      effective number of independent mechanisms.
%   C. Actual trial/kernel-integral distributions and IDR/IDQ/IDS designs.
%   D. Five 50-ms epochs to compare integrate-then-max (primary) with
%      max-then-integrate.  The prestep period remains excluded.
%
% Primary outputs
% ---------------
% results.rawTable      one row per stream, n, and operating mode
% results.summaryTable  symmetric directions averaged by circular offset
% results.operatingTable performance and step coherence by n/mode
% results.diagnostic    close-to-max binned counterfactual results
%
% When Config.SaveOutputs is true, MAT, CSV, and PDF outputs are written to
% Config.OutputDir.  All configuration fields and the calibrated fixed step
% are stored in results for reproducibility.
%
% Usage
% -----
%   results = runMaxPoolingNSweep;
%   results = runMaxPoolingNSweep('NTrials', 100000, ...
%       'InternalNoiseSD', 1.5, 'Visible', 'on');
%   cfg = maxPoolingDefaults; cfg.NValues = 2:10;
%   results = runMaxPoolingNSweep(cfg);
%
% See also RUNPAIREDPROBEOFFSETSWEEP, MAXPOOLINGDEFAULTS,
% SIMULATESPARSEMAXPOOLING, CALIBRATEMAXPOOLINGSTEP, FITPROBITIRLS.

cfg = maxPoolingDefaults();
cfg = parseConfig(cfg, varargin{:});

if ~exist(cfg.OutputDir, 'dir') && cfg.SaveOutputs
    mkdir(cfg.OutputDir);
end

% Calibrate the fixed step once at the reference n using independent trials.
refSeed = cfg.Seed + 100000 + cfg.ReferenceN;
refBank = makeSparseTrialBank(cfg.ReferenceN, cfg.NCalibrationTrials, cfg, refSeed);
fixedStep = calibrateMaxPoolingStep(refBank, cfg.TargetPerformance, cfg);

rawPieces = cell(2*numel(cfg.NValues), 1);
summaryPieces = cell(2*numel(cfg.NValues), 1);
operatingPieces = cell(2*numel(cfg.NValues), 1);
diagnosticSim = [];
piece = 0;

for ii = 1:numel(cfg.NValues)
    n = cfg.NValues(ii);
    calibrationBank = makeSparseTrialBank(n, cfg.NCalibrationTrials, cfg, ...
        cfg.Seed + 200000 + n);
    matchedStep = calibrateMaxPoolingStep(calibrationBank, ...
        cfg.TargetPerformance, cfg);

    analysisBank = makeSparseTrialBank(n, cfg.NTrials, cfg, ...
        cfg.Seed + 300000 + n);
    modeNames = {'Fixed', 'Matched'};
    modeSteps = [fixedStep, matchedStep];

    for mm = 1:2
        piece = piece + 1;
        sim = simulateSparseMaxPooling(analysisBank, modeSteps(mm), cfg);
        [rawPieces{piece}, summaryPieces{piece}, operatingPieces{piece}] = ...
            tablesFromSimulation(sim, modeNames{mm});

        if n == cfg.DiagnosticN && strcmp(modeNames{mm}, cfg.DiagnosticMode)
            diagnosticSim = sim;
        end
    end
end

results = struct();
results.config = cfg;
results.fixedStepCoherence = fixedStep;
results.rawTable = vertcat(rawPieces{1:piece});
results.summaryTable = vertcat(summaryPieces{1:piece});
results.operatingTable = vertcat(operatingPieces{1:piece});

if isempty(diagnosticSim)
    warning('Diagnostic n/mode was not simulated; using the first matched condition.');
    n = cfg.NValues(1);
    bank = makeSparseTrialBank(n, cfg.NTrials, cfg, cfg.Seed + 400000 + n);
    step = calibrateMaxPoolingStep(bank, cfg.TargetPerformance, cfg);
    diagnosticSim = simulateSparseMaxPooling(bank, step, cfg);
end
results.diagnostic = computeCloseToMaxDiagnostic(diagnosticSim, cfg);

results.figureHandles = plotSparseNSweepResults(results);

if cfg.SaveOutputs
    writetable(results.rawTable, fullfile(cfg.OutputDir, 'SparseNSweep_RawStreams.csv'));
    writetable(results.summaryTable, fullfile(cfg.OutputDir, 'SparseNSweep_OffsetSummary.csv'));
    writetable(results.operatingTable, fullfile(cfg.OutputDir, 'SparseNSweep_OperatingPoints.csv'));
    returnedFigureHandles = results.figureHandles;
    results.figureHandles = []; % Do not serialize live graphics objects.
    save(fullfile(cfg.OutputDir, 'SparseNSweep_Results.mat'), 'results', '-v7.3');
    results.figureHandles = returnedFigureHandles;
    saveSimulationFigures(results.figureHandles, cfg.OutputDir);
end
end

function cfg = parseConfig(cfg, varargin)
if numel(varargin) == 1 && isstruct(varargin{1})
    supplied = varargin{1};
    names = fieldnames(supplied);
    for k = 1:numel(names)
        if ~isfield(cfg, names{k})
            error('Unknown configuration field: %s', names{k});
        end
        cfg.(names{k}) = supplied.(names{k});
    end
elseif mod(numel(varargin), 2) == 0
    for k = 1:2:numel(varargin)
        name = char(varargin{k});
        if ~isfield(cfg, name)
            error('Unknown configuration field: %s', name);
        end
        cfg.(name) = varargin{k+1};
    end
else
    error('Use a configuration struct or name/value pairs.');
end

validateattributes(cfg.NValues, {'numeric'}, {'vector','integer','>=',2});
validateattributes(cfg.NTrials, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.NCalibrationTrials, {'numeric'}, {'scalar','integer','positive'});
validateattributes(cfg.TargetPerformance, {'numeric'}, {'scalar','>',0.5,'<',1});
validateattributes(cfg.TuningSigmaDeg, {'numeric'}, {'scalar','positive'});
validateattributes(cfg.ExternalNoiseSD, {'numeric'}, {'scalar','nonnegative'});
validateattributes(cfg.InternalNoiseSD, {'numeric'}, {'scalar','positive'});
validateattributes(cfg.StreamPedestal, {'numeric'}, {'vector','finite'});
end

function [rawTable, summaryTable, operatingTable] = tablesFromSimulation(sim, modeName)
n = sim.n;
modeCol = repmat(string(modeName), n, 1);
rawTable = table(repmat(n,n,1), modeCol, (1:n)', sim.directionsDeg(:), ...
    sim.offsetDeg(:), sim.fit.beta(:), sim.fit.se(:), ...
    sim.fit.beta(:)./sim.fit.beta(1), sim.winProbabilityChanged(:), ...
    repmat(sim.performance,n,1), repmat(sim.stepCoherence,n,1), ...
    'VariableNames', {'N','Mode','StreamIndex','DirectionDeg','OffsetDeg', ...
    'Gain','GainSE','NormalizedGain','AlignedMechanismWinProbability', ...
    'Performance','StepCoherence'});

offsets = unique(sim.offsetDeg(:));
nGroups = numel(offsets);
gain = nan(nGroups,1); gainSE = gain; normalized = gain; normalizedSE = gain;
winProb = gain; nStreams = gain;
covGain = sim.fit.covariance;
b0 = sim.fit.beta(1);
varB0 = covGain(1,1);

for g = 1:nGroups
    idx = find(abs(sim.offsetDeg(:)-offsets(g)) < 1e-9);
    c = zeros(n,1); c(idx) = 1/numel(idx);
    gain(g) = c' * sim.fit.beta;
    varG = max(c' * covGain * c, 0);
    gainSE(g) = sqrt(varG);
    covGB0 = c' * covGain(:,1);
    normalized(g) = gain(g)/b0;
    varRatio = varG/b0^2 + gain(g)^2*varB0/b0^4 ...
        - 2*gain(g)*covGB0/b0^3;
    normalizedSE(g) = sqrt(max(varRatio,0));
    winProb(g) = mean(sim.winProbabilityChanged(idx));
    nStreams(g) = numel(idx);
end

summaryTable = table(repmat(n,nGroups,1), repmat(string(modeName),nGroups,1), ...
    offsets, nStreams, gain, gainSE, normalized, normalizedSE, winProb, ...
    repmat(sim.performance,nGroups,1), repmat(sim.stepCoherence,nGroups,1), ...
    'VariableNames', {'N','Mode','OffsetDeg','NumSymmetricStreams','Gain', ...
    'GainSE','NormalizedGain','NormalizedGainSE', ...
    'AlignedMechanismWinProbability','Performance','StepCoherence'});

operatingTable = table(n, string(modeName), sim.performance, sim.stepCoherence, ...
    sim.fit.intercept, sim.fit.logLikelihood, sim.fit.converged, ...
    'VariableNames', {'N','Mode','Performance','StepCoherence', ...
    'ProbitIntercept','LogLikelihood','ProbitConverged'});
end
