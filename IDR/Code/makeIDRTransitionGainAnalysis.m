function analysis = makeIDRTransitionGainAnalysis(varargin)
% makeIDRTransitionGainAnalysis  Test transient noise responses in IDR choices.
%
% This analysis asks whether short responses to actual coherence transitions
% explain choice beyond the established temporally weighted sustained-noise
% predictors.  It is a strict extension of the INC/change-side gain model:
%
%   baseline:
%     dCeff = max(dC + gPref(parent)*IPref + gProbe*IProbe, 0)
%
%   transition model:
%     dCeff = max(dC + gPref(parent)*IPref + gProbe*IProbe + ...
%                         hPref*SPref + hProbe*SProbe, 0)
%
% I is the existing leave-one-parent-session-out weighted sustained
% coherence.  S is a signed transition waveform: every actual displayed
% coherence change contributes its change in coherence, followed by an
% exponential decay.  S is integrated with the same temporal weights as I.
%
% Important conventions
% ---------------------
% * Only certified INC trials with valid outcomes and coherence noise enter.
% * Only the change-side patch is analyzed, matching the current gain fit.
% * Stored event times are mapped to video frames with floor(), matching the
%   acquisition-era IDR renderer and extractPatchNoiseMatrices.
% * Events mapping to the same video frame are collapsed to their final
%   value before transitions are calculated.  The first displayed frame is
%   an initial state, not a transition.
% * Probe coherence is multiplied by two for paired yoked probes below 180
%   deg, and by one for the single 180-deg probe, matching the existing
%   effective-probe predictor.
% * Signal-step insertions do not count as noise transitions because the
%   stored noise values do not change.
% * The primary transient time constant is fixed a priori at 25 ms.
%   Sensitivity analyses use 12.5 and 50 ms by default.
% * Parent-session bootstrap confidence intervals are computed only for the
%   primary 25-ms model; cross-validation is computed for every time constant.
% * Cross-validation folds are assigned by parent recording session, never
%   by trial.
%
% Outputs
% -------
% Data/AcrossOffsetSummaries/GainAnalysis/
%   IDRTransitionGainAnalysis_<Animal>.mat
%   IDRTransitionGainAnalysis_<Animal>_allTau.csv
%   IDRTransitionGainAnalysis_<Animal>_primary.csv
%
% Plots/AcrossProbes/GainAnalysis/
%   IDRTransitionGainAnalysis_<Animal>.pdf
%
% Examples
% --------
%   A = makeIDRTransitionGainAnalysis('Animal','Meetz');
%   A = makeIDRTransitionGainAnalysis('Animal','Meetz','Replace',true);
%   A = makeIDRTransitionGainAnalysis('Animal','Meetz','NBootstrap',0);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(p, 'TauMS', [12.5 25 50], ...
  @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && all(x > 0));
addParameter(p, 'PrimaryTauMS', 25, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'CVFolds', 10, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2);
addParameter(p, 'CVSeed', 1729, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'NBootstrap', 200, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'BootstrapSeed', 2718, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'GainLimit', 5, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Bin180Into179', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'MinTrialsPerOffset', 100, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

tauMS = unique(double(opts.TauMS(:)'));
if ~any(abs(tauMS - opts.PrimaryTauMS) < 1e-10)
  error('makeIDRTransitionGainAnalysis:PrimaryTauMissing', ...
    'PrimaryTauMS (%g) must be included in TauMS.', opts.PrimaryTauMS);
end

root = domainFolder(mfilename('fullpath'));
dataRoot = fullfile(root, 'Data');
gainFolder = validFolder(fullfile(dataRoot, 'AcrossOffsetSummaries', 'GainAnalysis'));
plotFolder = validFolder(fullfile(root, 'Plots', 'AcrossProbes', 'GainAnalysis'));
animalTag = char(string(opts.Animal));
stem = sprintf('IDRTransitionGainAnalysis_%s', animalTag);
matPath = fullfile(gainFolder, [stem '.mat']);
allTauCSVPath = fullfile(gainFolder, [stem '_allTau.csv']);
primaryCSVPath = fullfile(gainFolder, [stem '_primary.csv']);
pdfPath = fullfile(plotFolder, [stem '.pdf']);

if isfile(matPath) && ~opts.Replace
  S = load(matPath, 'analysis');
  analysis = S.analysis;
  if opts.Verbose
    fprintf('IDR transition analysis already exists; loaded %s\n', matPath);
  end
  return;
end

weightPath = fullfile(dataRoot, 'AcrossOffsetSummaries', ...
  sprintf('BetaWeights_%s.mat', animalTag));
prefFitPath = fullfile(gainFolder, ...
  sprintf('IDRGainAcrossSessionPrefFit_%s.mat', animalTag));
if ~isfile(weightPath)
  error('makeIDRTransitionGainAnalysis:MissingWeights', ...
    'Missing BetaWeights file: %s', weightPath);
end
if ~isfile(prefFitPath)
  error('makeIDRTransitionGainAnalysis:MissingPreferredFit', ...
    'Missing preferred-gain fit: %s', prefFitPath);
end
W = load(weightPath, 'weightData');
P = load(prefFitPath, 'fitData');
weightData = W.weightData;
prefFit = P.fitData;

sessionFolder = fullfile(dataRoot, 'FullSessions');
[selectedFiles, ~] = selectAnalysisFiles(sessionFolder, 'Animal', opts.Animal);
if isempty(selectedFiles)
  error('makeIDRTransitionGainAnalysis:NoSessions', ...
    'No full IDR sessions passed selectAnalysisFiles.');
end

if opts.Verbose
  fprintf('Building sustained and transition predictors from %d candidate sessions...\n', ...
    numel(selectedFiles));
end
[trialTable, sessionTable] = buildTrialTable(selectedFiles, weightData, ...
  prefFit, tauMS, opts);
if isempty(trialTable)
  error('makeIDRTransitionGainAnalysis:NoTrials', ...
    'No eligible INC coherence-noise trials were found.');
end

offsets = unique(trialTable.probeDirDeg);
offsets = sort(offsets(:));
rows = repmat(emptyResultRow(), numel(offsets) * numel(tauMS), 1);
fits = cell(numel(offsets), numel(tauMS));
rowIndex = 0;

for iOffset = 1:numel(offsets)
  offset = offsets(iOffset);
  offsetMask = trialTable.probeDirDeg == offset;
  Dall = tableToData(trialTable(offsetMask,:));
  if numel(Dall.correct) < opts.MinTrialsPerOffset
    if opts.Verbose
      fprintf('Skipping offset %g: only %d eligible trials.\n', ...
        offset, numel(Dall.correct));
    end
    continue;
  end

  baseline = fitBaseline(Dall, opts.GainLimit);
  baselineCV = crossValidate(Dall, 'baseline', opts, NaN);

  for iTau = 1:numel(tauMS)
    rowIndex = rowIndex + 1;
    tau = tauMS(iTau);
    D = selectTau(Dall, iTau);
    augmented = fitAugmented(D, opts.GainLimit);
    augmentedCV = crossValidate(D, 'augmented', opts, tau);
    if abs(tau - opts.PrimaryTauMS) < 1e-10
      bootstrap = bootstrapAugmented(D, opts);
    else
      bootstrap = emptyBootstrap(opts.NBootstrap);
    end

    row = emptyResultRow();
    row.probeDirDeg = offset;
    row.tauMS = tau;
    row.isPrimaryTau = abs(tau - opts.PrimaryTauMS) < 1e-10;
    row.nTrials = numel(D.correct);
    row.nSessions = numel(unique(D.sessionIndex));
    row.nCorrect = sum(D.correct);
    row.fractionCorrect = mean(D.correct);
    row.gProbeBaseline = baseline.params(1);
    row.nllBaseline = baseline.nll;
    row.gProbeAugmented = augmented.params(1);
    row.hPref = augmented.params(2);
    row.hProbe = augmented.params(3);
    row.nllAugmented = augmented.nll;
    row.deltaNLLInSample = baseline.nll - augmented.nll;
    row.deltaGProbe = augmented.params(1) - baseline.params(1);
    row.fractionGProbeRemaining = safeRatio(augmented.params(1), baseline.params(1));
    row.cvNLLBaseline = baselineCV.nll;
    row.cvNLLAugmented = augmentedCV.nll;
    row.cvDeltaNLL = baselineCV.nll - augmentedCV.nll;
    row.cvDeltaNLLPerTrial = row.cvDeltaNLL / row.nTrials;
    row.corrIPrefSPref = corrFinite(D.iPref, D.sPref);
    row.corrIProbeSProbe = corrFinite(D.iProbe, D.sProbe);
    row.partialCorrIProbeSProbe = partialCorrFinite(D.iProbe, D.sProbe, ...
      [D.iPref, D.sPref]);
    row.hPrefCILow = bootstrap.hPrefCI(1);
    row.hPrefCIHigh = bootstrap.hPrefCI(2);
    row.hProbeCILow = bootstrap.hProbeCI(1);
    row.hProbeCIHigh = bootstrap.hProbeCI(2);
    row.gProbeAugCILow = bootstrap.gProbeCI(1);
    row.gProbeAugCIHigh = bootstrap.gProbeCI(2);
    row.bootstrapSuccessful = bootstrap.nSuccessful;
    rows(rowIndex) = row;

    fits{iOffset,iTau} = struct( ...
      'probeDirDeg', offset, ...
      'tauMS', tau, ...
      'baseline', baseline, ...
      'augmented', augmented, ...
      'baselineCV', baselineCV, ...
      'augmentedCV', augmentedCV, ...
      'bootstrap', bootstrap);

    if opts.Verbose
      fprintf(['  %5g deg, tau %5g ms: n=%5d, gP %.3f -> %.3f, ' ...
        'hP %.3f, CV DeltaNLL %.3f\n'], offset, tau, row.nTrials, ...
        row.gProbeBaseline, row.gProbeAugmented, row.hProbe, row.cvDeltaNLL);
    end
  end
end

rows = rows(1:rowIndex);
resultTable = struct2table(rows);
resultTable = sortrows(resultTable, {'probeDirDeg','tauMS'});
primaryTable = resultTable(resultTable.isPrimaryTau,:);

metadata = struct();
metadata.version = 1;
metadata.analysisName = 'IDRTransitionGainAnalysis';
metadata.primaryHypothesis = ['Short signed responses to actual displayed coherence transitions ' ...
  'explain IDR choices beyond the established sustained-noise predictors.'];
metadata.primaryTauMS = opts.PrimaryTauMS;
metadata.tauMS = tauMS;
metadata.transientKernel = 'causal exponential; q(k)=exp(-frameMS/tau)*q(k-1)+DeltaC(k)';
metadata.transitionDefinition = ['actual change in final displayed frame coherence after collapsing ' ...
  'all stored events that map to the same acquisition-era video frame'];
metadata.initialFrameConvention = 'initial displayed state is not counted as a transition';
metadata.probeConvention = ['effective yoked-pair coherence: 2x single probe stream for offsets ' ...
  'below 180 deg; 1x for 180 deg'];
metadata.sustainedWeighting = ['parent-session leave-one-out preferred-noise weights from ' ...
  'BetaWeights, normalized to sum to one over the step'];
metadata.crossValidation = sprintf('%d folds assigned by parent recording session', opts.CVFolds);
metadata.bootstrap = sprintf(['%d parent-session cluster resamples within offset; ' ...
  'primary tau only'], opts.NBootstrap);
metadata.weightPath = weightPath;
metadata.prefFitPath = prefFitPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

analysis = struct();
analysis.version = 1;
analysis.metadata = metadata;
analysis.sessionTable = sessionTable;
analysis.trialTable = trialTable;
analysis.resultTable = resultTable;
analysis.primaryTable = primaryTable;
analysis.fits = fits;

save(matPath, 'analysis', '-v7.3');
writetable(resultTable, allTauCSVPath);
writetable(primaryTable, primaryCSVPath);
makeSummaryFigure(primaryTable, pdfPath, animalTag, opts.PrimaryTauMS);

if opts.Verbose
  fprintf('Saved MAT: %s\n', matPath);
  fprintf('Saved primary summary: %s\n', primaryCSVPath);
  fprintf('Saved plot: %s\n', pdfPath);
end
end

% =========================================================================
function [T, STable] = buildTrialTable(selectedFiles, W, PF, tauMS, opts)

allTables = cell(numel(selectedFiles),1);
sRows = repmat(emptySessionRow(), numel(selectedFiles), 1);
nUsed = 0;

weightNames = baseNames(W.sessionFileNames);
fitNames = baseNames(PF.sessionTable.fileName);

for iFile = 1:numel(selectedFiles)
  filePath = selectedFiles{iFile};
  [~, sessionName] = fileparts(filePath);
  sessionBase = lower(string(sessionName));
  weightRow = find(weightNames == sessionBase);
  fitRow = find(fitNames == sessionBase);

  sRow = emptySessionRow();
  sRow.fileName = string(sessionName) + ".mat";
  sRow.filePath = string(filePath);
  sRow.weightRow = scalarOrNaN(weightRow);
  sRow.prefFitRow = scalarOrNaN(fitRow);

  if numel(weightRow) ~= 1 || numel(fitRow) ~= 1
    sRow.included = false;
    sRow.message = sprintf('Weight matches=%d, preferred-fit matches=%d.', ...
      numel(weightRow), numel(fitRow));
    sRows(iFile) = sRow;
    continue;
  end

  L = load(filePath, 'trials', 'sessionHeader');
  if ~isfield(L, 'trials') || ~isfield(L, 'sessionHeader')
    sRow.included = false;
    sRow.message = "Missing trials or sessionHeader.";
    sRows(iFile) = sRow;
    continue;
  end

  weights = double(W.leaveOneOutWeights(weightRow,:));
  stepTMS = double(W.stepTMS(:)');
  if abs(sum(weights) - 1) > 1e-12
    error('makeIDRTransitionGainAnalysis:BadWeights', ...
      'Weights do not sum to one for %s.', sessionName);
  end

  alpha = double(PF.sessionTable.alphaSession(fitRow));
  gPref = double(PF.sessionTable.gPrefSession(fitRow));
  beta = double(PF.params.betaShared);
  lapse = double(PF.params.lapseShared);

  H = L.sessionHeader;
  frameRateHz = headerScalar(H, 'frameRateHz');
  preStepMS = headerScalar(H, 'preStepMS');
  stepMS = headerScalar(H, 'stepMS');
  frameMS = 1000 / frameRateHz;
  nFrames = round((preStepMS + stepMS) / frameMS);
  tMS = (0:nFrames-1) * frameMS;
  stepMask = tMS >= preStepMS & tMS < preStepMS + stepMS;
  if numel(stepTMS) ~= sum(stepMask) || ...
      max(abs(stepTMS - tMS(stepMask))) > 1e-6
    error('makeIDRTransitionGainAnalysis:WeightTimeMismatch', ...
      'BetaWeights step times do not match %s.', sessionName);
  end

  trialRows = repmat(emptyTrialRow(numel(tauMS)), numel(L.trials), 1);
  nTrialRows = 0;
  for iTrial = 1:numel(L.trials)
    tr = L.trials{iTrial};
    if ~eligibleTrial(tr)
      continue;
    end
    D = tr.trial.data;
    probeDir = double(D.probeDirDeg);
    if probeDir <= 0 || probeDir > 180
      continue;
    end
    outputOffset = probeDir;
    if probeDir == 180 && opts.Bin180Into179
      outputOffset = 179;
    end
    nProbeStreams = 1 + double(probeDir < 180);

    times = numericData(tr.changeTimesMS);
    prefValues = numericData(tr.changePrefCohsPC);
    probeValues = nProbeStreams .* numericData(tr.changeProbeCohsPC);
    [prefFrames, probeFrames] = displayedFrameSequences(times, prefValues, ...
      probeValues, nFrames, frameMS);

    iPref = sum(weights .* prefFrames(stepMask));
    iProbe = sum(weights .* probeFrames(stepMask));
    sPref = nan(1,numel(tauMS));
    sProbe = nan(1,numel(tauMS));
    for iTau = 1:numel(tauMS)
      qPref = transitionWaveform(prefFrames, frameMS, tauMS(iTau));
      qProbe = transitionWaveform(probeFrames, frameMS, tauMS(iTau));
      sPref(iTau) = sum(weights .* qPref(stepMask));
      sProbe(iTau) = sum(weights .* qProbe(stepMask));
    end

    nTrialRows = nTrialRows + 1;
    r = emptyTrialRow(numel(tauMS));
    r.sessionIndex = iFile;
    r.sessionName = string(sessionName);
    r.parentTrialIndex = iTrial;
    r.probeDirDeg = outputOffset;
    r.originalProbeDirDeg = probeDir;
    r.nYokedProbeStreams = nProbeStreams;
    r.correct = double(numericData(tr.trialEnd) == 0);
    r.deltaC = double(D.stepCohPC) - double(D.baseCohPC);
    r.alpha = alpha;
    r.gPref = gPref;
    r.beta = beta;
    r.lapse = lapse;
    r.iPref = iPref;
    r.iProbe = iProbe;
    r.sPref = sPref;
    r.sProbe = sProbe;
    r.prefStart = prefFrames(find(stepMask,1,'first'));
    r.prefEnd = prefFrames(find(stepMask,1,'last'));
    r.probeStart = probeFrames(find(stepMask,1,'first'));
    r.probeEnd = probeFrames(find(stepMask,1,'last'));
    r.nPrefTransitionsStep = countTransitions(prefFrames, stepMask);
    r.nProbeTransitionsStep = countTransitions(probeFrames, stepMask);
    trialRows(nTrialRows) = r;
  end

  trialRows = trialRows(1:nTrialRows);
  if nTrialRows > 0
    nUsed = nUsed + 1;
    allTables{nUsed} = struct2table(trialRows);
    sRow.included = true;
    sRow.nTrials = nTrialRows;
    sRow.message = "ok";
  else
    sRow.message = "No eligible trials.";
  end
  sRows(iFile) = sRow;
end

allTables = allTables(1:nUsed);
if isempty(allTables)
  T = table();
else
  T = vertcat(allTables{:});
end
STable = struct2table(sRows);
end

% =========================================================================
function tf = eligibleTrial(tr)
tf = false;
required = {'trial','trialEnd','trialCertify','changeTimesMS', ...
  'changePrefCohsPC','changeProbeCohsPC'};
if ~all(isfield(tr, required)) || ~isfield(tr.trial, 'data')
  return;
end
D = tr.trial.data;
needed = {'changeIndex','cohNoise','probeDirDeg','baseCohPC','stepCohPC'};
if ~all(isfield(D, needed))
  return;
end
outcome = numericData(tr.trialEnd);
certify = numericData(tr.trialCertify);
tf = isscalar(outcome) && ismember(outcome,[0 1]) && ...
  isscalar(certify) && certify == 0 && ...
  double(D.changeIndex) == 1 && logical(D.cohNoise);
end

% =========================================================================
function [prefFrames, probeFrames] = displayedFrameSequences(times, pref, probe, nFrames, frameMS)
times = double(times(:));
pref = double(pref(:));
probe = double(probe(:));
n = min([numel(times), numel(pref), numel(probe)]);
if n == 0
  error('makeIDRTransitionGainAnalysis:EmptySequence', ...
    'Eligible noise trial has an empty change-side sequence.');
end
times = times(1:n);
pref = pref(1:n);
probe = probe(1:n);
if any(diff(times) < 0) || times(1) ~= 0
  error('makeIDRTransitionGainAnalysis:BadSequenceTimes', ...
    'Noise event times must begin at zero and be nondecreasing.');
end

% Acquisition-era mapping: integer event time is floored to a video-frame
% index.  Assign sequentially so the final event in a frame supersedes any
% earlier event in that same frame.
eventFrame = floor(times ./ frameMS) + 1;
eventFrame = max(eventFrame, 1);
keep = eventFrame <= nFrames;
eventFrame = eventFrame(keep);
pref = pref(keep);
probe = probe(keep);

prefAtFrame = nan(1,nFrames);
probeAtFrame = nan(1,nFrames);
for i = 1:numel(eventFrame)
  prefAtFrame(eventFrame(i)) = pref(i);
  probeAtFrame(eventFrame(i)) = probe(i);
end
eventFrames = find(isfinite(prefAtFrame) & isfinite(probeAtFrame));
if isempty(eventFrames) || eventFrames(1) ~= 1
  error('makeIDRTransitionGainAnalysis:NoInitialState', ...
    'Noise sequence does not define the first displayed video frame.');
end

prefFrames = nan(1,nFrames);
probeFrames = nan(1,nFrames);
for i = 1:numel(eventFrames)
  f0 = eventFrames(i);
  if i < numel(eventFrames)
    f1 = eventFrames(i+1)-1;
  else
    f1 = nFrames;
  end
  prefFrames(f0:f1) = prefAtFrame(f0);
  probeFrames(f0:f1) = probeAtFrame(f0);
end
if any(~isfinite(prefFrames)) || any(~isfinite(probeFrames))
  error('makeIDRTransitionGainAnalysis:IncompleteSequence', ...
    'Failed to reconstruct every displayed noise frame.');
end
end

% =========================================================================
function q = transitionWaveform(coherenceFrames, frameMS, tauMS)
decay = exp(-frameMS / tauMS);
delta = [0, diff(coherenceFrames)];
q = zeros(size(coherenceFrames));
for k = 2:numel(q)
  q(k) = decay * q(k-1) + delta(k);
end
end

% =========================================================================
function n = countTransitions(frames, stepMask)
first = find(stepMask,1,'first');
last = find(stepMask,1,'last');
% Include a transition at the first step frame relative to the preceding
% frame, but do not include a transition after the final step frame.
idx = max(first-1,1):last;
n = sum(diff(frames(idx)) ~= 0);
end

% =========================================================================
function D = tableToData(T)
D = struct();
vars = {'correct','deltaC','alpha','gPref','beta','lapse','iPref','iProbe', ...
  'sessionIndex','prefStart','prefEnd','probeStart','probeEnd', ...
  'nPrefTransitionsStep','nProbeTransitionsStep'};
for i = 1:numel(vars)
  D.(vars{i}) = double(T.(vars{i})(:));
end
D.sPrefAll = double(T.sPref);
D.sProbeAll = double(T.sProbe);
end

% =========================================================================
function D = selectTau(Dall, iTau)
D = Dall;
D.sPref = Dall.sPrefAll(:,iTau);
D.sProbe = Dall.sProbeAll(:,iTau);
end

% =========================================================================
function fit = fitBaseline(D, gainLimit)
theta0 = gainToTheta(0, gainLimit);
objective = @(theta) modelNLL(theta, D, 'baseline', gainLimit);
[theta, nll, exitflag, output] = runFminsearch(objective, theta0);
fit.params = gainLimit .* tanh(theta(:)');
fit.paramNames = {'gProbe'};
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% =========================================================================
function fit = fitAugmented(D, gainLimit)
theta0 = gainToTheta([0 0 0], gainLimit);
objective = @(theta) modelNLL(theta, D, 'augmented', gainLimit);
[theta, nll, exitflag, output] = runFminsearch(objective, theta0);
fit.params = gainLimit .* tanh(theta(:)');
fit.paramNames = {'gProbe','hPref','hProbe'};
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% =========================================================================
function [theta, nll, exitflag, output] = runFminsearch(objective, theta0)
options = optimset('Display','off','MaxIter',2000,'MaxFunEvals',10000, ...
  'TolX',1e-8,'TolFun',1e-8);
[theta, nll, exitflag, output] = fminsearch(objective, theta0, options);
end

% =========================================================================
function nll = modelNLL(theta, D, model, gainLimit)
params = gainLimit .* tanh(theta(:)');
p = modelProbability(params, D, model);
nll = bernoulliNLL(p, D.correct);
end

% =========================================================================
function p = modelProbability(params, D, model)
effective = D.deltaC + D.gPref .* D.iPref + params(1) .* D.iProbe;
if strcmp(model, 'augmented')
  effective = effective + params(2) .* D.sPref + params(3) .* D.sProbe;
end
effective = max(effective, 0);
p = 0.5 + (0.5 - D.lapse) .* ...
  (1 - exp(-((effective ./ D.alpha) .^ D.beta)));
p = min(max(p,1e-12),1-1e-12);
end

% =========================================================================
function nll = bernoulliNLL(p, correct)
nll = -sum(correct .* log(p) + (1-correct) .* log(1-p));
if ~isfinite(nll), nll = realmax; end
end

% =========================================================================
function cv = crossValidate(D, model, opts, tau)
sessions = unique(D.sessionIndex);
nFolds = min(round(opts.CVFolds), numel(sessions));
cv = struct('nll',NaN,'foldNLL',nan(nFolds,1),'foldNTrials',zeros(nFolds,1), ...
  'tauMS',tau,'model',model);
if nFolds < 2
  return;
end

oldRNG = rng;
cleanup = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
rng(opts.CVSeed,'twister');
sessions = sessions(randperm(numel(sessions)));
foldBySession = mod((0:numel(sessions)-1)',nFolds)+1;

for iFold = 1:nFolds
  testSessions = sessions(foldBySession == iFold);
  test = ismember(D.sessionIndex,testSessions);
  train = ~test;
  if ~any(test) || ~any(train)
    continue;
  end
  Dtrain = subsetData(D,train);
  Dtest = subsetData(D,test);
  if strcmp(model,'baseline')
    fit = fitBaseline(Dtrain,opts.GainLimit);
  else
    fit = fitAugmented(Dtrain,opts.GainLimit);
  end
  p = modelProbability(fit.params,Dtest,model);
  cv.foldNLL(iFold) = bernoulliNLL(p,Dtest.correct);
  cv.foldNTrials(iFold) = sum(test);
end
cv.nll = sum(cv.foldNLL,'omitnan');
end

% =========================================================================
function B = bootstrapAugmented(D, opts)
B = emptyBootstrap(opts.NBootstrap);
if opts.NBootstrap <= 0
  return;
end
sessions = unique(D.sessionIndex);
if numel(sessions) < 2
  return;
end
oldRNG = rng;
cleanup = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
rng(opts.BootstrapSeed,'twister');

for iBoot = 1:opts.NBootstrap
  sampled = sessions(randi(numel(sessions),numel(sessions),1));
  indices = cell(numel(sampled),1);
  for k = 1:numel(sampled)
    indices{k} = find(D.sessionIndex == sampled(k));
  end
  idx = vertcat(indices{:});
  Db = subsetData(D,idx);
  try
    fit = fitAugmented(Db,opts.GainLimit);
    if fit.exitflag > 0 && all(isfinite(fit.params))
      B.params(iBoot,:) = fit.params;
    end
  catch
    % Retain NaNs for failed resamples and report nSuccessful.
  end
end
ok = all(isfinite(B.params),2);
B.nSuccessful = sum(ok);
if B.nSuccessful >= max(20,ceil(0.5*opts.NBootstrap))
  B.gProbeCI = percentileLocal(B.params(ok,1),[2.5 97.5]);
  B.hPrefCI = percentileLocal(B.params(ok,2),[2.5 97.5]);
  B.hProbeCI = percentileLocal(B.params(ok,3),[2.5 97.5]);
end
end

% =========================================================================
function B = emptyBootstrap(nRequested)
B = struct('nRequested',nRequested,'nSuccessful',0, ...
  'params',nan(nRequested,3),'gProbeCI',[NaN NaN], ...
  'hPrefCI',[NaN NaN],'hProbeCI',[NaN NaN]);
end

% =========================================================================
function out = subsetData(D, index)
out = D;
fields = fieldnames(D);
n = numel(D.correct);
for i = 1:numel(fields)
  x = D.(fields{i});
  if (isnumeric(x) || islogical(x)) && size(x,1) == n
    out.(fields{i}) = x(index,:);
  end
end
end

% =========================================================================
function makeSummaryFigure(T, pdfPath, animalTag, tauMS)
if isempty(T), return; end
x = T.probeDirDeg;
fig = figure('Color','w','Position',[100 100 1100 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(x,T.gProbeBaseline,'ko-','LineWidth',1.3,'MarkerFaceColor','k'); hold on;
errorbar(x,T.gProbeAugmented, ...
  T.gProbeAugmented-T.gProbeAugCILow, ...
  T.gProbeAugCIHigh-T.gProbeAugmented, ...
  'o-','Color',[0.1 0.45 0.8],'MarkerFaceColor',[0.1 0.45 0.8], ...
  'LineWidth',1.3);
yline(0,'k:'); xlabel('Probe offset (deg)'); ylabel('Sustained probe gain');
legend({'Baseline','With transitions'},'Location','best','Box','off'); box off;

nexttile;
errorbar(x,T.hProbe,T.hProbe-T.hProbeCILow,T.hProbeCIHigh-T.hProbe, ...
  'o-','Color',[0.75 0.2 0.15],'MarkerFaceColor',[0.75 0.2 0.15], ...
  'LineWidth',1.3); hold on;
yline(0,'k:'); xlabel('Probe offset (deg)'); ylabel('Probe transition gain h_P'); box off;

nexttile;
bar(x,T.cvDeltaNLLPerTrial,'FaceColor',[0.35 0.65 0.4]); hold on;
yline(0,'k:'); xlabel('Probe offset (deg)'); ylabel('CV Delta NLL / trial'); box off;

nexttile;
plot(x,T.corrIProbeSProbe,'o-','LineWidth',1.3); hold on;
plot(x,T.partialCorrIProbeSProbe,'s-','LineWidth',1.3);
yline(0,'k:'); xlabel('Probe offset (deg)'); ylabel('Predictor correlation');
legend({'corr(I_P,S_P)','partial corr'},'Location','best','Box','off'); box off;

sgtitle(sprintf('IDR transition-gain analysis: %s, tau = %g ms',animalTag,tauMS));
exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end

% =========================================================================
function r = emptyResultRow()
r = struct('probeDirDeg',NaN,'tauMS',NaN,'isPrimaryTau',false, ...
  'nTrials',NaN,'nSessions',NaN,'nCorrect',NaN,'fractionCorrect',NaN, ...
  'gProbeBaseline',NaN,'nllBaseline',NaN,'gProbeAugmented',NaN, ...
  'hPref',NaN,'hProbe',NaN,'nllAugmented',NaN,'deltaNLLInSample',NaN, ...
  'deltaGProbe',NaN,'fractionGProbeRemaining',NaN, ...
  'cvNLLBaseline',NaN,'cvNLLAugmented',NaN,'cvDeltaNLL',NaN, ...
  'cvDeltaNLLPerTrial',NaN,'corrIPrefSPref',NaN,'corrIProbeSProbe',NaN, ...
  'partialCorrIProbeSProbe',NaN,'hPrefCILow',NaN,'hPrefCIHigh',NaN, ...
  'hProbeCILow',NaN,'hProbeCIHigh',NaN,'gProbeAugCILow',NaN, ...
  'gProbeAugCIHigh',NaN,'bootstrapSuccessful',NaN);
end

% =========================================================================
function r = emptySessionRow()
r = struct('fileName',"",'filePath',"",'weightRow',NaN,'prefFitRow',NaN, ...
  'included',false,'nTrials',0,'message',"");
end

% =========================================================================
function r = emptyTrialRow(nTau)
r = struct('sessionIndex',NaN,'sessionName',"",'parentTrialIndex',NaN, ...
  'probeDirDeg',NaN,'originalProbeDirDeg',NaN,'nYokedProbeStreams',NaN, ...
  'correct',NaN,'deltaC',NaN,'alpha',NaN,'gPref',NaN,'beta',NaN,'lapse',NaN, ...
  'iPref',NaN,'iProbe',NaN,'sPref',nan(1,nTau),'sProbe',nan(1,nTau), ...
  'prefStart',NaN,'prefEnd',NaN,'probeStart',NaN,'probeEnd',NaN, ...
  'nPrefTransitionsStep',NaN,'nProbeTransitionsStep',NaN);
end

% =========================================================================
function names = baseNames(paths)
paths = string(paths(:));
names = strings(size(paths));
for i = 1:numel(paths)
  [~,b] = fileparts(paths(i));
  names(i) = lower(string(b));
end
end

% =========================================================================
function x = numericData(x)
while isstruct(x) && isfield(x,'data')
  x = x.data;
end
x = double(x);
end

% =========================================================================
function value = headerScalar(H,fieldName)
if ~isfield(H,fieldName)
  error('makeIDRTransitionGainAnalysis:MissingHeaderField', ...
    'sessionHeader.%s is required.',fieldName);
end
value = numericData(H.(fieldName));
if ~isscalar(value) || ~isfinite(value)
  error('makeIDRTransitionGainAnalysis:BadHeaderField', ...
    'sessionHeader.%s must be a finite scalar.',fieldName);
end
end

% =========================================================================
function theta = gainToTheta(g,gainLimit)
z = max(min(g ./ gainLimit,0.95),-0.95);
theta = atanh(z);
end

% =========================================================================
function x = scalarOrNaN(x)
if numel(x) == 1, x = double(x); else, x = NaN; end
end

% =========================================================================
function r = safeRatio(a,b)
if isfinite(a) && isfinite(b) && abs(b) > 1e-12
  r = a/b;
else
  r = NaN;
end
end

% =========================================================================
function r = corrFinite(x,y)
ok = isfinite(x) & isfinite(y);
if sum(ok) < 3 || std(x(ok)) == 0 || std(y(ok)) == 0
  r = NaN;
else
  C = corrcoef(x(ok),y(ok));
  r = C(1,2);
end
end

% =========================================================================
function r = partialCorrFinite(x,y,Z)
ok = isfinite(x) & isfinite(y) & all(isfinite(Z),2);
if sum(ok) < size(Z,2)+3
  r = NaN;
  return;
end
X = [ones(sum(ok),1),Z(ok,:)];
rx = x(ok)-X*(X\x(ok));
ry = y(ok)-X*(X\y(ok));
r = corrFinite(rx,ry);
end

% =========================================================================
function q = percentileLocal(x,pct)
x = sort(x(isfinite(x)));
q = nan(size(pct));
if isempty(x), return; end
for i = 1:numel(pct)
  pos = 1 + (numel(x)-1)*pct(i)/100;
  lo = floor(pos); hi = ceil(pos);
  if lo == hi
    q(i) = x(lo);
  else
    q(i) = x(lo)+(pos-lo)*(x(hi)-x(lo));
  end
end
end
