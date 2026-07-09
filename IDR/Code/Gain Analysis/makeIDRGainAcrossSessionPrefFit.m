function fitData = makeIDRGainAcrossSessionPrefFit(varargin)
% makeIDRGainAcrossSessionPrefFit_v1
% Joint preferred-noise gain fit for IDR Inc/ChangeSide trials.
%
% Reads existing full-session increment products in:
%   Data/FullSessions/BetaAnalysis/*.mat
%
% Uses the existing leave-one-session-out preferred-noise kernel weights in:
%   Data/AcrossOffsetSummaries/BetaWeights_<Animal>.mat
%
% Fits the model requested for the stabilized preferred-noise stage:
%
%   deltaC_eff = max(deltaC + gPref_session * nPref, 0)
%
%   P(correct) = 0.5 + (0.5 - lapseShared) * ...
%      (1 - exp(-(deltaC_eff / alpha_session)^betaShared))
%
% Parameters:
%   params.alphaBySession
%   params.gPrefBySession
%   params.betaShared
%   params.lapseShared
%
% alpha and nPref are in percent coherence. gPref is dimensionless.
%
% This first version uses alternating optimization by default, which is much
% better conditioned than applying fminsearch to all ~2*nSessions+2
% parameters at once.
%
% Outputs:
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefFit_<Animal>.mat
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainAcrossSessionPrefFit_<Animal>_sessionTable.csv

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'InitialBeta', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'InitialLapse', 0.02, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 0.5);
addParameter(p, 'MaxLapse', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 0.5);
addParameter(p, 'GainLimit', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'MinTrials', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'MaxOuterIter', 25, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'TolNLL', 1e-6, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'DoFullPolish', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

if opts.InitialLapse > opts.MaxLapse
  error('makeIDRGainAcrossSessionPrefFit:BadInitialLapse', ...
    'InitialLapse must be <= MaxLapse.');
end

root = domainFolder(mfilename('fullpath'));
sessionFolder = fullfile(root, 'Data', 'FullSessions', 'BetaAnalysis');
outputFolder = validFolder(fullfile(root, 'Data', 'AcrossOffsetSummaries', 'GainAnalysis'));

animalTag = char(string(opts.Animal));
matPath = fullfile(outputFolder, sprintf('IDRGainAcrossSessionPrefFit_%s.mat', animalTag));
% csvPath = fullfile(outputFolder, sprintf('IDRGainAcrossSessionPrefFit_%s_sessionTable.csv', animalTag));

if isfile(matPath) && ~opts.Replace
  if opts.Verbose
    fprintf('IDR across-session preferred-gain fit already exists; loading %s\n', matPath);
  end
  S = load(matPath, 'fitData');
  fitData = S.fitData;
  return;
end

[selectedFiles, fileInfo] = selectAnalysisFiles(sessionFolder, 'Animal', opts.Animal);
if isempty(selectedFiles)
  error('makeIDRGainAcrossSessionPrefFit:NoSelectedFiles', ...
    'No BetaAnalysis session files passed selectAnalysisFiles.');
end

weightPath = fullfile(root, 'Data', 'AcrossOffsetSummaries', sprintf('BetaWeights_%s.mat', animalTag));
if ~isfile(weightPath) && ~strcmpi(animalTag, 'All')
  allWeightPath = fullfile(root, 'Data', 'AcrossOffsetSummaries', 'BetaWeights_All.mat');
  if isfile(allWeightPath)
    weightPath = allWeightPath;
  end
end
if ~isfile(weightPath)
  error('makeIDRGainAcrossSessionPrefFit:MissingWeights', ...
    'Could not find BetaWeights_<Animal>.mat or BetaWeights_All.mat. Run makeBetaKernel first.');
end
W = load(weightPath, 'weightData');
weightData = W.weightData;

rawRows = repmat(emptyRawRow(), numel(selectedFiles), 1);
trialSessionIndex = [];
deltaC = [];
nPref = [];
correct = [];

for iFile = 1:numel(selectedFiles)
  filePath = selectedFiles{iFile};
  if opts.Verbose
    fprintf('Loading preferred-gain trials: %s\n', fileInfo.fileName{iFile});
  end

  row = emptyRawRow();
  row.animal = string(fileInfo.animal{iFile});
  row.sessionName = string(stripExtension(fileInfo.fileName{iFile}));
  row.fileName = string(fileInfo.fileName{iFile});
  row.filePath = string(filePath);

  try
    S = load(filePath, 'sessionNoise');
    N = S.sessionNoise;
    % requireSessionNoiseFields(N, filePath);

    c = double(N.signalCohPC(:));
    y = double(N.trialOutcome(:) == 0);
    [np, prefInfo] = computePreferredNoisePredictor(N, weightData, filePath);
    np = double(np(:));

    use = isfinite(c) & c >= 0 & isfinite(y) & isfinite(np);
    c = c(use);
    y = y(use);
    np = np(use);

    row.nTrials = double(getWithDefault(N, 'nTrials', numel(N.trialOutcome)));
    row.nValidTrials = numel(c);
    row.nIncChangeSideTrials = row.nValidTrials;
    row.nPreferredNoiseTrials = sum(logical(N.hasPreferredNoise(:)) & use(:));
    row.nNoPreferredNoiseTrials = sum(~logical(N.hasPreferredNoise(:)) & use(:));
    row.parentWeightRow = prefInfo.parentWeightRow;
    row.leaveOneOutKernelSum = prefInfo.leaveOneOutKernelSum;
    row.prefNoiseMean = mean(np);
    row.prefNoiseSD = std(np);
    row.prefNoiseMin = min(np);
    row.prefNoiseMax = max(np);
    row.fractionCorrect = mean(y);
    row.minCoherence = min(c);
    row.maxCoherence = max(c);
    row.nCoherenceLevels = numel(unique(c));
    row.fitIncluded = row.nValidTrials >= opts.MinTrials && row.nCoherenceLevels >= 2;
    row.message = "ok";

    if row.fitIncluded
      fitIndex = height(struct2table(rawRows(1:iFile))) ; %#ok<NASGU>
      sessionNumber = iFile;
      trialSessionIndex = [trialSessionIndex; repmat(sessionNumber, numel(c), 1)]; %#ok<AGROW>
      deltaC = [deltaC; c]; %#ok<AGROW>
      nPref = [nPref; np]; %#ok<AGROW>
      correct = [correct; y]; %#ok<AGROW>
    else
      row.message = "Too few trials or coherence levels for stabilized fit.";
    end
  catch ME
    row.fitIncluded = false;
    row.message = string(ME.message);
  end
  rawRows(iFile) = row;
end

rawSessionTable = struct2table(rawRows);
includedRows = find(rawSessionTable.fitIncluded);
nSessions = numel(includedRows);
if nSessions < 1
  error('makeIDRGainAcrossSessionPrefFit:NoUsableSessions', ...
    'No sessions were usable for the across-session preferred-gain fit.');
end

% Remap trialSessionIndex from original file order to contiguous included-session order.
oldToNew = nan(numel(selectedFiles), 1);
oldToNew(includedRows) = 1:nSessions;
trialSessionIndex = oldToNew(trialSessionIndex);
if any(~isfinite(trialSessionIndex))
  error('makeIDRGainAcrossSessionPrefFit:SessionIndexBug', ...
    'Internal session-index remapping failed.');
end
trialSessionIndex = double(trialSessionIndex(:));

% Starting values.
betaShared = opts.InitialBeta;
lapseShared = opts.InitialLapse;
alphaBySession = nan(nSessions, 1);
gPrefBySession = ones(nSessions, 1);

for iSession = 1:nSessions
  m = trialSessionIndex == iSession;
  alphaBySession(iSession) = initialAlpha(deltaC(m));
end

history = table();
prevNLL = Inf;
for iOuter = 1:opts.MaxOuterIter
  % Step 1: fit each session's alpha and gPref with beta/lapse fixed.
  sessionExitflag = nan(nSessions, 1);
  sessionMessage = strings(nSessions, 1);
  for iSession = 1:nSessions
    m = trialSessionIndex == iSession;
    f = fitSessionAlphaGain(deltaC(m), nPref(m), correct(m), betaShared, lapseShared, alphaBySession(iSession), ...
      gPrefBySession(iSession), opts.GainLimit);
    alphaBySession(iSession) = f.alpha;
    gPrefBySession(iSession) = f.gPref;
    sessionExitflag(iSession) = f.exitflag;
    sessionMessage(iSession) = string(f.message);
  end

  % Step 2: fit shared beta and lapse with session alpha/gPref fixed.
  sharedFit = fitSharedBetaLapse(deltaC, nPref, correct, trialSessionIndex, ...
    alphaBySession, gPrefBySession, betaShared, lapseShared, opts.MaxLapse);
  betaShared = sharedFit.betaShared;
  lapseShared = sharedFit.lapseShared;

  totalNLL = preferredAcrossNLL(deltaC, nPref, correct, trialSessionIndex, ...
    alphaBySession, gPrefBySession, betaShared, lapseShared);

  hrow = table(iOuter, totalNLL, betaShared, lapseShared, ...
    min(alphaBySession), max(alphaBySession), mean(gPrefBySession), std(gPrefBySession), sharedFit.exitflag, ...
    string(sharedFit.message), 'VariableNames', {'iteration','nll','betaShared','lapseShared', ...
    'minAlpha','maxAlpha','meanGPref','sdGPref','sharedExitflag','sharedMessage'});
  history = [history; hrow]; %#ok<AGROW>

  if opts.Verbose
    fprintf('  iter %02d: NLL %.6f, beta %.4g, lapse %.4g, mean gPref %.4g\n', ...
      iOuter, totalNLL, betaShared, lapseShared, mean(gPrefBySession));
  end

  if abs(prevNLL - totalNLL) < opts.TolNLL
    break;
  end
  prevNLL = totalNLL;
end

fullPolish = struct('used', opts.DoFullPolish, 'exitflag', NaN, 'message', '', 'nll', NaN);
if opts.DoFullPolish
  polished = fullPolishAllParams(deltaC, nPref, correct, trialSessionIndex, ...
    alphaBySession, gPrefBySession, betaShared, lapseShared, opts.MaxLapse, opts.GainLimit);
  alphaBySession = polished.alphaBySession;
  gPrefBySession = polished.gPrefBySession;
  betaShared = polished.betaShared;
  lapseShared = polished.lapseShared;
  fullPolish.exitflag = polished.exitflag;
  fullPolish.message = polished.message;
  fullPolish.nll = polished.nll;
end

% Final likelihood and per-session contributions.
totalNLL = preferredAcrossNLL(deltaC, nPref, correct, trialSessionIndex, ...
  alphaBySession, gPrefBySession, betaShared, lapseShared);

sessionTable = rawSessionTable(includedRows, :);
sessionTable.fitSessionIndex = (1:nSessions)';
sessionTable.alphaSession = alphaBySession;
sessionTable.gPrefSession = gPrefBySession;
sessionTable.betaShared = repmat(betaShared, nSessions, 1);
sessionTable.lapseShared = repmat(lapseShared, nSessions, 1);
sessionTable.nllContribution = nan(nSessions, 1);
sessionTable.meanPredictedPCorrect = nan(nSessions, 1);

pHat = preferredAcrossPcorrect(deltaC, nPref, trialSessionIndex, ...
  alphaBySession, gPrefBySession, betaShared, lapseShared);
for iSession = 1:nSessions
  m = trialSessionIndex == iSession;
  sessionTable.nllContribution(iSession) = bernoulliNLL(pHat(m), correct(m));
  sessionTable.meanPredictedPCorrect(iSession) = mean(pHat(m));
end

trialTable = table();
trialTable.sessionIndex = trialSessionIndex;
trialTable.sessionName = sessionTable.sessionName(trialSessionIndex);
trialTable.deltaC = deltaC;
trialTable.nPref = nPref;
trialTable.correct = logical(correct);
trialTable.pHat = pHat;
trialTable.deltaCEff = max(deltaC + gPrefBySession(trialSessionIndex) .* nPref, 0);

params = struct();
params.alphaBySession = alphaBySession;
params.gPrefBySession = gPrefBySession;
params.betaShared = betaShared;
params.lapseShared = lapseShared;

metadata = struct();
metadata.version = 1;
metadata.analysisName = 'IDRGainAcrossSessionPrefFit';
metadata.method = ['Alternating stabilized preferred-noise gain fit; ' ...
  'alpha and preferred gain are session-specific; beta and lapse are shared.'];
metadata.model = ['deltaC_eff = max(deltaC + gPref_session*nPref, 0); ' ...
  'Pcorrect = 0.5 + (0.5-lapseShared)*(1-exp(-(deltaC_eff/alpha_session)^betaShared))'];
metadata.deltaCUnits = 'percent coherence';
metadata.noisePredictorUnits = 'percent coherence';
metadata.gPrefUnits = 'dimensionless coherence-equivalent gain';
metadata.sourceFolder = sessionFolder;
metadata.weightPath = weightPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

fitData = struct();
fitData.version = 1;
fitData.params = params;
fitData.nll = totalNLL;
fitData.nTrials = height(trialTable);
fitData.nSessions = nSessions;
fitData.sessionTable = sessionTable;
fitData.rawSessionTable = rawSessionTable;
fitData.trialTable = trialTable;
fitData.history = history;
fitData.fullPolish = fullPolish;
fitData.metadata = metadata;

save(matPath, 'fitData', '-v7.3');
% writetable(sessionTable, csvPath);

if opts.Verbose
  fprintf('Saved MAT: %s\n', matPath);
  % fprintf('Saved CSV: %s\n', csvPath);
end
end

% -------------------------------------------------------------------------
function row = emptyRawRow()
row = struct();
row.animal = "";
row.sessionName = "";
row.fileName = "";
row.filePath = "";
row.nTrials = NaN;
row.nValidTrials = NaN;
row.nIncChangeSideTrials = NaN;
row.nPreferredNoiseTrials = NaN;
row.nNoPreferredNoiseTrials = NaN;
row.nCoherenceLevels = NaN;
row.minCoherence = NaN;
row.maxCoherence = NaN;
row.fractionCorrect = NaN;
row.prefNoiseMean = NaN;
row.prefNoiseSD = NaN;
row.prefNoiseMin = NaN;
row.prefNoiseMax = NaN;
row.parentWeightRow = NaN;
row.leaveOneOutKernelSum = NaN;
row.fitIncluded = false;
row.message = "";
end

% -------------------------------------------------------------------------
% function requireSessionNoiseFields(N, filePath)
% required = {'trialOutcome', 'signalCohPC', 'sessionHeader', ...
%   'hasPreferredNoise', 'noiseTimesMS', 'noiseCohsPC'};
% missing = required(~isfield(N, required));
% if ~isempty(missing)
%   error('makeIDRGainAcrossSessionPrefFit:MissingFields', ...
%     '%s is missing sessionNoise fields: %s', filePath, strjoin(missing, ', '));
% end
% end

% -------------------------------------------------------------------------
function [nPref, info] = computePreferredNoisePredictor(N, W, sessionPath)
weights = preferredWeightsForSession(W, sessionPath);
H = N.sessionHeader;
frameRateHz = headerScalarLocal(H, 'frameRateHz');
preStepMS = headerScalarLocal(H, 'preStepMS');
stepMS = headerScalarLocal(H, 'stepMS');
msPerFrame = 1000 / frameRateHz;
nFrames = round((preStepMS + stepMS) / msPerFrame);
tMS = (0:nFrames-1) * msPerFrame;
stepMask = tMS >= preStepMS & tMS < preStepMS + stepMS;
stepTMS = tMS(stepMask);

if numel(weights.weights) ~= numel(stepTMS)
  error('makeIDRGainAcrossSessionPrefFit:WeightLengthMismatch', ...
    'Weight count (%d) does not match step-frame count (%d).', numel(weights.weights), numel(stepTMS));
end
if max(abs(stepTMS(:) - weights.stepTMS(:))) > 1e-6
  error('makeIDRGainAcrossSessionPrefFit:WeightTimeMismatch', ...
    'Session step-frame times do not match BetaWeights.stepTMS.');
end
if abs(sum(weights.weights) - 1) > 1e-12
  error('makeIDRGainAcrossSessionPrefFit:BadWeightNormalization', ...
    'Preferred-noise weights do not sum to one.');
end

nTrials = numel(N.trialOutcome);
nPref = zeros(nTrials, 1);
hasPref = logical(N.hasPreferredNoise(:));

for iTrial = 1:nTrials
  if ~hasPref(iTrial)
    nPref(iTrial) = 0;
    continue;
  end
  times = double(N.noiseTimesMS{iTrial}(:)');
  values = double(N.noiseCohsPC{iTrial}(:)');
  sampled = samplePiecewiseConstant(times, values, stepTMS);
  nPref(iTrial) = sum(sampled(:)' .* weights.weights(:)');
end

info = struct();
info.nPreferredNoiseTrials = sum(hasPref);
info.nNoPreferredNoiseTrials = sum(~hasPref);
info.parentWeightRow = weights.parentWeightRow;
info.leaveOneOutKernelSum = weights.leaveOneOutKernelSum;
end

% -------------------------------------------------------------------------
function out = preferredWeightsForSession(W, sessionPath)
required = {'sessionFileNames', 'leaveOneOutWeights', 'stepTMS', 'leaveOneOutKernelSum'};
missing = required(~isfield(W, required));
if ~isempty(missing)
  error('makeIDRGainAcrossSessionPrefFit:MissingWeightFields', ...
    'weightData is missing fields: %s', strjoin(missing, ', '));
end
[~, sessionBase, sessionExt] = fileparts(sessionPath);
sessionName = [sessionBase sessionExt];
weightNames = string(W.sessionFileNames(:));
weightBases = strings(size(weightNames));
for i = 1:numel(weightNames)
  [~, b, e] = fileparts(weightNames(i));
  weightBases(i) = b + e;
end
row = find(strcmpi(weightBases, sessionName));
if numel(row) ~= 1
  error('makeIDRGainAcrossSessionPrefFit:WeightSessionMatch', ...
    'Could not uniquely match %s to BetaWeights sessionFileNames (%d matches).', sessionName, numel(row));
end
out = struct();
out.parentWeightRow = row;
out.stepTMS = double(W.stepTMS(:));
out.weights = double(W.leaveOneOutWeights(row, :));
out.leaveOneOutKernelSum = double(W.leaveOneOutKernelSum(row));
end

% -------------------------------------------------------------------------
function sampled = samplePiecewiseConstant(timesMS, valuesPC, queryTimesMS)
if isempty(timesMS) || isempty(valuesPC)
  error('makeIDRGainAcrossSessionPrefFit:EmptyNoiseSequence', ...
    'Preferred-noise trials must contain nonempty times and values.');
end
if numel(timesMS) ~= numel(valuesPC)
  error('makeIDRGainAcrossSessionPrefFit:NoiseLengthMismatch', ...
    'Noise times and values must have equal length for preferred-noise trials.');
end
sampled = nan(size(queryTimesMS));
for i = 1:numel(queryTimesMS)
  idx = find(timesMS <= queryTimesMS(i), 1, 'last');
  if isempty(idx)
    error('makeIDRGainAcrossSessionPrefFit:NoPrecedingNoiseSample', ...
      'No piecewise-constant noise sample precedes query time %.6g.', queryTimesMS(i));
  end
  sampled(i) = valuesPC(idx);
end
end

% -------------------------------------------------------------------------
function alpha0 = initialAlpha(deltaC)
positiveC = deltaC(deltaC > 0 & isfinite(deltaC));
if isempty(positiveC)
  alpha0 = 25;
else
  alpha0 = median(positiveC);
end
if ~isfinite(alpha0) || alpha0 <= 0
  alpha0 = max(positiveC) / 2;
end
if ~isfinite(alpha0) || alpha0 <= 0
  alpha0 = 25;
end
end

% -------------------------------------------------------------------------
function fit = fitSessionAlphaGain(deltaC, nPref, correct, betaShared, lapseShared, alpha0, g0, gainLimit)
theta0 = [log(max(alpha0, eps)), gainToTheta(g0, gainLimit)];
objective = @(theta) sessionAlphaGainNLL(theta, deltaC, nPref, correct, betaShared, lapseShared, gainLimit);
options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000, ...
  'TolX', 1e-7, 'TolFun', 1e-7);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, options);
fit.alpha = exp(thetaHat(1));
fit.gPref = gainLimit * tanh(thetaHat(2));
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function nll = sessionAlphaGainNLL(theta, deltaC, nPref, correct, betaShared, lapseShared, gainLimit)
alpha = exp(theta(1));
gPref = gainLimit * tanh(theta(2));
deltaCEff = max(deltaC + gPref .* nPref, 0);
p = weibullPcorrect(deltaCEff, alpha, betaShared, lapseShared);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function fit = fitSharedBetaLapse(deltaC, nPref, correct, sessionIndex, alphaBySession, gPrefBySession, beta0, lapse0, maxLapse)
theta0 = [log(max(beta0, eps)), lapseToTheta(lapse0, maxLapse)];
objective = @(theta) sharedBetaLapseNLL(theta, deltaC, nPref, correct, ...
  sessionIndex, alphaBySession, gPrefBySession, maxLapse);
options = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 3000, ...
  'TolX', 1e-8, 'TolFun', 1e-8);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, options);
fit.betaShared = exp(thetaHat(1));
fit.lapseShared = maxLapse ./ (1 + exp(-thetaHat(2)));
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function nll = sharedBetaLapseNLL(theta, deltaC, nPref, correct, sessionIndex, alphaBySession, gPrefBySession, maxLapse)
betaShared = exp(theta(1));
lapseShared = maxLapse ./ (1 + exp(-theta(2)));
p = preferredAcrossPcorrect(deltaC, nPref, sessionIndex, alphaBySession, ...
  gPrefBySession, betaShared, lapseShared);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function polished = fullPolishAllParams(deltaC, nPref, correct, sessionIndex, alphaBySession, gPrefBySession, betaShared, lapseShared, maxLapse, gainLimit)
theta0 = packAllTheta(alphaBySession, gPrefBySession, betaShared, lapseShared, maxLapse, gainLimit);
objective = @(theta) allParamNLL(theta, deltaC, nPref, correct, sessionIndex, numel(alphaBySession), maxLapse, gainLimit);
options = optimset('Display', 'off', 'MaxIter', 5000, 'MaxFunEvals', 20000, ...
  'TolX', 1e-7, 'TolFun', 1e-7);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, options);
[alphaBySession, gPrefBySession, betaShared, lapseShared] = unpackAllTheta(thetaHat, numel(alphaBySession), maxLapse, gainLimit);
polished = struct();
polished.alphaBySession = alphaBySession;
polished.gPrefBySession = gPrefBySession;
polished.betaShared = betaShared;
polished.lapseShared = lapseShared;
polished.nll = nll;
polished.exitflag = exitflag;
polished.message = output.message;
end

% -------------------------------------------------------------------------
function theta = packAllTheta(alphaBySession, gPrefBySession, betaShared, lapseShared, maxLapse, gainLimit)
theta = [log(alphaBySession(:)); arrayfun(@(g) gainToTheta(g, gainLimit), gPrefBySession(:)); ...
  log(betaShared); lapseToTheta(lapseShared, maxLapse)];
end

% -------------------------------------------------------------------------
function [alphaBySession, gPrefBySession, betaShared, lapseShared] = unpackAllTheta(theta, nSessions, maxLapse, gainLimit)
alphaBySession = exp(theta(1:nSessions));
gTheta = theta(nSessions+1:2*nSessions);
gPrefBySession = gainLimit .* tanh(gTheta);
betaShared = exp(theta(2*nSessions+1));
lapseShared = maxLapse ./ (1 + exp(-theta(2*nSessions+2)));
end

% -------------------------------------------------------------------------
function nll = allParamNLL(theta, deltaC, nPref, correct, sessionIndex, nSessions, maxLapse, gainLimit)
[alphaBySession, gPrefBySession, betaShared, lapseShared] = unpackAllTheta(theta, nSessions, maxLapse, gainLimit);
nll = preferredAcrossNLL(deltaC, nPref, correct, sessionIndex, ...
  alphaBySession, gPrefBySession, betaShared, lapseShared);
end

% -------------------------------------------------------------------------
function nll = preferredAcrossNLL(deltaC, nPref, correct, sessionIndex, alphaBySession, gPrefBySession, betaShared, lapseShared)
p = preferredAcrossPcorrect(deltaC, nPref, sessionIndex, alphaBySession, ...
  gPrefBySession, betaShared, lapseShared);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function p = preferredAcrossPcorrect(deltaC, nPref, sessionIndex, alphaBySession, gPrefBySession, betaShared, lapseShared)
alpha = alphaBySession(sessionIndex);
gPref = gPrefBySession(sessionIndex);
deltaCEff = max(deltaC + gPref .* nPref, 0);
p = weibullPcorrect(deltaCEff, alpha, betaShared, lapseShared);
end

% -------------------------------------------------------------------------
function p = weibullPcorrect(deltaC, alpha, beta, lapse)
p = 0.5 + (0.5 - lapse) .* (1 - exp(-((deltaC ./ alpha) .^ beta)));
end

% -------------------------------------------------------------------------
function nll = bernoulliNLL(p, correct)
epsP = 1e-12;
p = min(max(p, epsP), 1 - epsP);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));
if ~isfinite(nll)
  nll = realmax;
end
end

% -------------------------------------------------------------------------
function theta = gainToTheta(g, gainLimit)
z = max(min(g ./ gainLimit, 0.95), -0.95);
theta = atanh(z);
end

% -------------------------------------------------------------------------
function theta = lapseToTheta(lapse, maxLapse)
if maxLapse == 0
  theta = -Inf;
  return;
end
z = max(min(lapse ./ maxLapse, 1 - 1e-6), 1e-6);
theta = log(z ./ (1 - z));
end

% -------------------------------------------------------------------------
function out = getWithDefault(S, fieldName, defaultValue)
if isfield(S, fieldName)
  out = S.(fieldName);
else
  out = defaultValue;
end
end

% -------------------------------------------------------------------------
function stem = stripExtension(fileName)
[~, stem] = fileparts(char(fileName));
end

% -------------------------------------------------------------------------
function value = headerScalarLocal(H, fieldName)
if ~isfield(H, fieldName)
  error('makeIDRGainAcrossSessionPrefFit:MissingHeaderField', ...
    'sessionHeader.%s is required.', fieldName);
end
value = H.(fieldName);
while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  error('makeIDRGainAcrossSessionPrefFit:BadHeaderScalar', ...
    'sessionHeader.%s must be a numeric scalar.', fieldName);
end
value = double(value);
end
