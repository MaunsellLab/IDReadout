function diagnosticTable = makeIDRGainSessionDiagnostics(varargin)
% makeIDRGainSessionDiagnostics  Per-session Weibull/gain diagnostics for IDR.
%
% Reads existing full-session increment products in:
%
%   Data/FullSessions/BetaAnalysis/*.mat
%
% and writes an across-session diagnostic table to:
%
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainSessionDiagnostics_<Animal>.mat
%   Data/AcrossOffsetSummaries/GainAnalysis/IDRGainSessionDiagnostics_<Animal>.csv
%
% Scope:
%   - IDR Inc / ChangeSide trials only, as represented by the existing
%     BetaAnalysis sessionNoise products.
%   - No no-change trials are included.
%   - Psychometric floor is fixed at 0.5 at zero increment coherence.
%   - The pure Weibull fits ignore all noise predictors.
%   - The preferred-gain fits use only the change-side preferred-noise
%     stream. Probe-noise streams are not used here.
%
% Pure Weibull model:
%   P(correct) = 0.5 + (0.5 - lambda) * ... 
%      (1 - exp(-(deltaC / alpha)^beta))
%
% Preferred-gain model:
%   deltaC_eff = max(deltaC + gPref * nPref, 0)
%   P(correct) = Weibull(deltaC_eff; alpha, beta, lambda)
%
% where deltaC and nPref are in percent coherence, so gPref is dimensionless.
%
% Usage:
%   T = makeIDRGainSessionDiagnostics();
%   T = makeIDRGainSessionDiagnostics('Animal','Meetz','Replace',true);
%   T = makeIDRGainSessionDiagnostics('FixedBeta',2,'MaxLapse',0.05);
%   T = makeIDRGainSessionDiagnostics('DoPreferredGain',false);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FixedBeta', 2, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'MaxLapse', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x < 0.5);
addParameter(p, 'MinTrials', 20, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p, 'DoPreferredGain', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'GainLimit', 5, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
sessionFolder = fullfile(root, 'Data', 'FullSessions', 'BetaAnalysis');
outputFolder = validFolder(fullfile(root, 'Data', 'AcrossOffsetSummaries', 'GainAnalysis'));

animalTag = char(string(opts.Animal));
matPath = fullfile(outputFolder, sprintf('IDRGainSessionDiagnostics_%s.mat', animalTag));
csvPath = fullfile(outputFolder, sprintf('IDRGainSessionDiagnostics_%s.csv', animalTag));

if isfile(matPath) && isfile(csvPath) && ~opts.Replace
  if opts.Verbose
    fprintf('IDR gain diagnostics already exist; loading %s\n', matPath);
  end
  S = load(matPath, 'diagnosticTable');
  diagnosticTable = S.diagnosticTable;
  return;
end

[selectedFiles, fileInfo] = selectAnalysisFiles(sessionFolder, 'Animal', opts.Animal);
if isempty(selectedFiles)
  error('makeIDRGainSessionDiagnostics:NoSelectedFiles', ...
    'No BetaAnalysis session files passed selectAnalysisFiles.');
end

weightData = [];
weightPath = '';
if opts.DoPreferredGain
  weightPath = fullfile(root, 'Data', 'AcrossOffsetSummaries', sprintf('BetaWeights_%s.mat', animalTag));
  if ~isfile(weightPath) && ~strcmpi(animalTag, 'All')
    % Many existing runs create the all-session weight file. Fall back to it
    % for convenience, but record the actual path in metadata.
    allWeightPath = fullfile(root, 'Data', 'AcrossOffsetSummaries', 'BetaWeights_All.mat');
    if isfile(allWeightPath)
      weightPath = allWeightPath;
    end
  end
  if ~isfile(weightPath)
    error('makeIDRGainSessionDiagnostics:MissingWeights', ...
      ['Preferred-gain fits require BetaWeights_<Animal>.mat in Data/AcrossOffsetSummaries. ' ...
       'Run makeBetaKernel first, or call with DoPreferredGain=false.']);
  end
  W = load(weightPath, 'weightData');
  weightData = W.weightData;
end

rows = repmat(emptyRow(), numel(selectedFiles), 1);

for iFile = 1:numel(selectedFiles)
  filePath = selectedFiles{iFile};
  if opts.Verbose
    fprintf('Fitting IDR gain Weibull diagnostic: %s\n', fileInfo.fileName{iFile});
  end

  row = emptyRow();
  row.fileName = string(fileInfo.fileName{iFile});
  row.filePath = string(filePath);
  row.sessionName = string(stripExtension(fileInfo.fileName{iFile}));
  row.animal = string(fileInfo.animal{iFile});

  try
    S = load(filePath, 'sessionHeader', 'sessionNoise');
    if ~isfield(S, 'sessionNoise')
      error('File does not contain sessionNoise.');
    end
    N = S.sessionNoise;
    requireSessionNoiseFields(N, filePath, opts.DoPreferredGain);

    correct = double(N.trialOutcome(:) == 0);
    deltaC = double(N.signalCohPC(:));
    valid = isfinite(deltaC) & deltaC >= 0 & isfinite(correct);

    row.nTrials = double(getWithDefault(N, 'nTrials', numel(deltaC)));
    row.nValidTrials = sum(valid);

    % Existing BetaAnalysis files are already restricted to valid increment
    % trials from the change-side patch. Keep this field explicit so later
    % products can detect if the selection contract changes.
    row.nIncChangeSideTrials = row.nValidTrials;

    deltaC = deltaC(valid);
    correct = correct(valid);

    [cohLevels, nByCoh, nCorrectByCoh, pCorrectByCoh] = summarizeByCoherence(deltaC, correct);
    row.coherenceLevels = string(mat2str(cohLevels(:)', 6));
    row.nByCoh = string(mat2str(nByCoh(:)'));
    row.nCorrectByCoh = string(mat2str(nCorrectByCoh(:)'));
    row.pCorrectByCoh = string(mat2str(pCorrectByCoh(:)', 4));
    row.nCoherenceLevels = numel(cohLevels);
    row.minCoherence = min(deltaC);
    row.maxCoherence = max(deltaC);
    row.fractionCorrect = mean(correct);

    if row.nValidTrials < opts.MinTrials || numel(cohLevels) < 2
      row.message = "Too few trials or coherence levels for diagnostic Weibull fit.";
      rows(iFile) = row;
      continue;
    end

    fitFree = fitWeibullDiagnostic(deltaC, correct, [], opts.MaxLapse);
    row.alphaFree = fitFree.alpha;
    row.betaFree = fitFree.beta;
    row.lapseFree = fitFree.lapse;
    row.nllFree = fitFree.nll;
    row.exitflagFree = fitFree.exitflag;
    row.messageFree = string(fitFree.message);

    fitFixed = fitWeibullDiagnostic(deltaC, correct, opts.FixedBeta, opts.MaxLapse);
    row.fixedBeta = opts.FixedBeta;
    row.alphaFixed = fitFixed.alpha;
    row.lapseFixed = fitFixed.lapse;
    row.nllFixed = fitFixed.nll;
    row.exitflagFixed = fitFixed.exitflag;
    row.messageFixed = string(fitFixed.message);

    if opts.DoPreferredGain
      [nPref, prefInfo] = computePreferredNoisePredictor(N, weightData, filePath);
      nPref = nPref(valid);

      row.nPreferredNoiseTrials = prefInfo.nPreferredNoiseTrials;
      row.nNoPreferredNoiseTrials = prefInfo.nNoPreferredNoiseTrials;
      row.nGainTrials = numel(nPref);
      row.prefNoiseMean = mean(nPref);
      row.prefNoiseSD = std(nPref);
      row.prefNoiseMin = min(nPref);
      row.prefNoiseMax = max(nPref);
      row.parentWeightRow = prefInfo.parentWeightRow;
      row.leaveOneOutKernelSum = prefInfo.leaveOneOutKernelSum;

      gainFree = fitPreferredGainDiagnostic(deltaC, nPref, correct, [], opts.MaxLapse, opts.GainLimit);
      row.gPrefFree = gainFree.gain;
      row.alphaGainFree = gainFree.alpha;
      row.betaGainFree = gainFree.beta;
      row.lapseGainFree = gainFree.lapse;
      row.nllGainFree = gainFree.nll;
      row.exitflagGainFree = gainFree.exitflag;
      row.messageGainFree = string(gainFree.message);

      gainFixed = fitPreferredGainDiagnostic(deltaC, nPref, correct, opts.FixedBeta, opts.MaxLapse, opts.GainLimit);
      row.gPrefFixed = gainFixed.gain;
      row.alphaGainFixed = gainFixed.alpha;
      row.lapseGainFixed = gainFixed.lapse;
      row.nllGainFixed = gainFixed.nll;
      row.exitflagGainFixed = gainFixed.exitflag;
      row.messageGainFixed = string(gainFixed.message);
    end

    row.fitUsable = row.exitflagFixed > 0 && isfinite(row.nllFixed);
    if opts.DoPreferredGain
      row.fitUsable = row.fitUsable && row.exitflagGainFixed > 0 && isfinite(row.nllGainFixed);
    end
    if row.fitUsable
      row.message = "ok";
    else
      row.message = "One or more diagnostic fits did not converge cleanly.";
    end
  catch ME
    row.fitUsable = false;
    row.message = string(ME.message);
  end

  rows(iFile) = row;
end

diagnosticTable = struct2table(rows);
metadata = struct();
metadata.version = 2;
metadata.analysisName = 'IDRGainSessionDiagnostics';
metadata.method = ['Per-session Inc/ChangeSide Weibull and preferred-noise gain diagnostics from existing ' ...
  'Data/FullSessions/BetaAnalysis sessionNoise products'];
metadata.model = ['Pcorrect = 0.5 + (0.5 - lambda) * ' ...
  '(1 - exp(-(deltaC_eff / alpha)^beta))'];
metadata.deltaCUnits = 'percent coherence';
metadata.noisePredictorUnits = 'percent coherence';
metadata.preferredGainModel = 'deltaC_eff = max(deltaC + gPref * nPref, 0)';
metadata.probeNoiseIncluded = false;
metadata.fixedBeta = opts.FixedBeta;
metadata.maxLapse = opts.MaxLapse;
metadata.gainLimit = opts.GainLimit;
metadata.sourceFolder = sessionFolder;
metadata.weightPath = weightPath;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

save(matPath, 'diagnosticTable', 'metadata', '-v7.3');
writetable(diagnosticTable, csvPath);

if opts.Verbose
  fprintf('Saved MAT: %s\n', matPath);
  fprintf('Saved CSV: %s\n', csvPath);
end
end

% -------------------------------------------------------------------------
function row = emptyRow()
row = struct();
row.animal = "";
row.sessionName = "";
row.fileName = "";
row.filePath = "";
row.nTrials = NaN;
row.nValidTrials = NaN;
row.nIncChangeSideTrials = NaN;
row.nCoherenceLevels = NaN;
row.coherenceLevels = "";
row.nByCoh = "";
row.nCorrectByCoh = "";
row.pCorrectByCoh = "";
row.minCoherence = NaN;
row.maxCoherence = NaN;
row.fractionCorrect = NaN;
row.alphaFree = NaN;
row.betaFree = NaN;
row.lapseFree = NaN;
row.nllFree = NaN;
row.exitflagFree = NaN;
row.messageFree = "";
row.fixedBeta = NaN;
row.alphaFixed = NaN;
row.lapseFixed = NaN;
row.nllFixed = NaN;
row.exitflagFixed = NaN;
row.messageFixed = "";
row.nPreferredNoiseTrials = NaN;
row.nNoPreferredNoiseTrials = NaN;
row.nGainTrials = NaN;
row.prefNoiseMean = NaN;
row.prefNoiseSD = NaN;
row.prefNoiseMin = NaN;
row.prefNoiseMax = NaN;
row.parentWeightRow = NaN;
row.leaveOneOutKernelSum = NaN;
row.gPrefFree = NaN;
row.alphaGainFree = NaN;
row.betaGainFree = NaN;
row.lapseGainFree = NaN;
row.nllGainFree = NaN;
row.exitflagGainFree = NaN;
row.messageGainFree = "";
row.gPrefFixed = NaN;
row.alphaGainFixed = NaN;
row.lapseGainFixed = NaN;
row.nllGainFixed = NaN;
row.exitflagGainFixed = NaN;
row.messageGainFixed = "";
row.fitUsable = false;
row.message = "";
end

% -------------------------------------------------------------------------
function requireSessionNoiseFields(N, filePath, doPreferredGain)
required = {'trialOutcome', 'signalCohPC'};
if doPreferredGain
  required = [required, {'sessionHeader', 'hasPreferredNoise', 'noiseTimesMS', 'noiseCohsPC'}];
end
missing = required(~isfield(N, required));
if ~isempty(missing)
  error('makeIDRGainSessionDiagnostics:MissingFields', ...
    '%s is missing sessionNoise fields: %s', filePath, strjoin(missing, ', '));
end
end

% -------------------------------------------------------------------------
function [levels, nByCoh, nCorrectByCoh, pCorrectByCoh] = summarizeByCoherence(deltaC, correct)
levels = unique(deltaC(:));
nByCoh = zeros(size(levels));
nCorrectByCoh = zeros(size(levels));
for i = 1:numel(levels)
  mask = deltaC == levels(i);
  nByCoh(i) = sum(mask);
  nCorrectByCoh(i) = sum(correct(mask));
end
pCorrectByCoh = nCorrectByCoh ./ nByCoh;
end

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
  error('makeIDRGainSessionDiagnostics:WeightLengthMismatch', ...
    'Weight count (%d) does not match step-frame count (%d).', numel(weights.weights), numel(stepTMS));
end
if max(abs(stepTMS(:) - weights.stepTMS(:))) > 1e-6
  error('makeIDRGainSessionDiagnostics:WeightTimeMismatch', ...
    'Session step-frame times do not match BetaWeights.stepTMS.');
end
if abs(sum(weights.weights) - 1) > 1e-12
  error('makeIDRGainSessionDiagnostics:BadWeightNormalization', ...
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
  error('makeIDRGainSessionDiagnostics:MissingWeightFields', ...
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
  error('makeIDRGainSessionDiagnostics:WeightSessionMatch', ...
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
  error('makeIDRGainSessionDiagnostics:EmptyNoiseSequence', ...
    'Preferred-noise trials must contain nonempty times and values.');
end
if numel(timesMS) ~= numel(valuesPC)
  error('makeIDRGainSessionDiagnostics:NoiseLengthMismatch', ...
    'Noise times and values must have equal length for preferred-noise trials.');
end
sampled = nan(size(queryTimesMS));
for i = 1:numel(queryTimesMS)
  idx = find(timesMS <= queryTimesMS(i), 1, 'last');
  if isempty(idx)
    error('makeIDRGainSessionDiagnostics:NoPrecedingNoiseSample', ...
      'No piecewise-constant noise sample precedes query time %.6g.', queryTimesMS(i));
  end
  sampled(i) = valuesPC(idx);
end
end

% -------------------------------------------------------------------------
function fit = fitWeibullDiagnostic(deltaC, correct, fixedBeta, maxLapse)
% Fit using fminsearch on unconstrained transformed parameters, avoiding a
% dependency on the Optimization Toolbox.

deltaC = double(deltaC(:));
correct = double(correct(:));

positiveC = deltaC(deltaC > 0 & isfinite(deltaC));
if isempty(positiveC)
  error('makeIDRGainSessionDiagnostics:NoPositiveCoherence', ...
    'At least one positive increment coherence is required.');
end

alpha0 = median(positiveC);
if ~isfinite(alpha0) || alpha0 <= 0
  alpha0 = max(positiveC) / 2;
end
if isempty(fixedBeta)
  theta0 = [log(alpha0), log(2), 0];
else
  theta0 = [log(alpha0), 0];
end

objective = @(theta) weibullNLL(theta, deltaC, correct, fixedBeta, maxLapse);
opts = optimset('Display', 'off', 'MaxIter', 5000, 'MaxFunEvals', 10000, ...
  'TolX', 1e-8, 'TolFun', 1e-8);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, opts);

[alpha, beta, lapse] = unpackTheta(thetaHat, fixedBeta, maxLapse);
fit = struct();
fit.alpha = alpha;
fit.beta = beta;
fit.lapse = lapse;
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function fit = fitPreferredGainDiagnostic(deltaC, nPref, correct, fixedBeta, maxLapse, gainLimit)
deltaC = double(deltaC(:));
nPref = double(nPref(:));
correct = double(correct(:));

positiveC = deltaC(deltaC > 0 & isfinite(deltaC));
if isempty(positiveC)
  error('makeIDRGainSessionDiagnostics:NoPositiveCoherence', ...
    'At least one positive increment coherence is required.');
end

alpha0 = median(positiveC);
if ~isfinite(alpha0) || alpha0 <= 0
  alpha0 = max(positiveC) / 2;
end
gain0 = 1;
gainTheta0 = atanh(max(min(gain0 / gainLimit, 0.95), -0.95));
if isempty(fixedBeta)
  theta0 = [log(alpha0), log(2), 0, gainTheta0];
else
  theta0 = [log(alpha0), 0, gainTheta0];
end

objective = @(theta) preferredGainNLL(theta, deltaC, nPref, correct, fixedBeta, maxLapse, gainLimit);
opts = optimset('Display', 'off', 'MaxIter', 10000, 'MaxFunEvals', 20000, ...
  'TolX', 1e-8, 'TolFun', 1e-8);
[thetaHat, nll, exitflag, output] = fminsearch(objective, theta0, opts);

[alpha, beta, lapse, gain] = unpackGainTheta(thetaHat, fixedBeta, maxLapse, gainLimit);
fit = struct();
fit.alpha = alpha;
fit.beta = beta;
fit.lapse = lapse;
fit.gain = gain;
fit.nll = nll;
fit.exitflag = exitflag;
fit.message = output.message;
end

% -------------------------------------------------------------------------
function nll = weibullNLL(theta, deltaC, correct, fixedBeta, maxLapse)
[alpha, beta, lapse] = unpackTheta(theta, fixedBeta, maxLapse);
p = weibullPcorrect(deltaC, alpha, beta, lapse);
nll = bernoulliNLL(p, correct);
end

% -------------------------------------------------------------------------
function nll = preferredGainNLL(theta, deltaC, nPref, correct, fixedBeta, maxLapse, gainLimit)
[alpha, beta, lapse, gain] = unpackGainTheta(theta, fixedBeta, maxLapse, gainLimit);
deltaCEff = max(deltaC + gain .* nPref, 0);
p = weibullPcorrect(deltaCEff, alpha, beta, lapse);
nll = bernoulliNLL(p, correct);
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
function [alpha, beta, lapse] = unpackTheta(theta, fixedBeta, maxLapse)
alpha = exp(theta(1));
if isempty(fixedBeta)
  beta = exp(theta(2));
  lapseTheta = theta(3);
else
  beta = fixedBeta;
  lapseTheta = theta(2);
end
lapse = maxLapse ./ (1 + exp(-lapseTheta));
end

% -------------------------------------------------------------------------
function [alpha, beta, lapse, gain] = unpackGainTheta(theta, fixedBeta, maxLapse, gainLimit)
alpha = exp(theta(1));
if isempty(fixedBeta)
  beta = exp(theta(2));
  lapseTheta = theta(3);
  gainTheta = theta(4);
else
  beta = fixedBeta;
  lapseTheta = theta(2);
  gainTheta = theta(3);
end
lapse = maxLapse ./ (1 + exp(-lapseTheta));
gain = gainLimit .* tanh(gainTheta);
end

% -------------------------------------------------------------------------
function p = weibullPcorrect(deltaC, alpha, beta, lapse)
p = 0.5 + (0.5 - lapse) .* (1 - exp(-((deltaC ./ alpha) .^ beta)));
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
  error('makeIDRGainSessionDiagnostics:MissingHeaderField', ...
    'sessionHeader.%s is required.', fieldName);
end
value = H.(fieldName);
while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  error('makeIDRGainSessionDiagnostics:BadHeaderScalar', ...
    'sessionHeader.%s must be a numeric scalar.', fieldName);
end
value = double(value);
end
