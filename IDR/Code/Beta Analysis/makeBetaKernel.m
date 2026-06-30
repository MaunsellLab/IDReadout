function kernelData = makeBetaKernel(varargin)
% makeBetaKernel  Construct pooled preferred-noise reverse-correlation kernel.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%
% Uses only trials with:
%   sessionNoise.hasPreferredNoise == true
%
% All valid increment trials remain in SessionData for later psychometric
% fitting, but no-preferred-noise trials do not enter the kernel.
%
% Usage:
%   kernelData = makeBetaKernel();
%   kernelData = makeBetaKernel(3);
%   kernelData = makeBetaKernel("IDReadout2_Meetz_20260610.mat");

p = inputParser;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
parse(p, varargin{:});
R = p.Results;

baseFolder = domainFolder(mfilename('fullpath'));
sessionDataFolder = char(fullfile(baseFolder, 'Data', 'FullSessions', 'BetaAnalysis'));
outputFolder = validFolder(fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries'));
plotFolder = validFolder(fullfile(baseFolder, 'Plots', 'AcrossProbes', 'Kernels'));
[selectedFiles, fileInfo] = selectAnalysisFiles({sessionDataFolder}, 'Animal', R.Animal);
nSessions = numel(selectedFiles);

initialized = false;
allNoise = [];
allOutcome = [];
allSessionIndex = [];

sessionKernelsPC = [];
sessionKernelSEMPC = [];
sessionMeanCorrectPC = [];
sessionMeanErrorPC = [];
sessionNCorrect = zeros(1, nSessions);
sessionNError = zeros(1, nSessions);
sessionNTrials = zeros(1, nSessions);
sessionNAllIncrementTrials = zeros(1, nSessions);
sessionSignalCohPC = cell(1, nSessions);
sessionPrefCohNoisePC = nan(1, nSessions);

for iSession = 1:nSessions
  load(selectedFiles{iSession}, 'sessionNoise');
  H = sessionNoise.sessionHeader;
  requiredFields = {'hasPreferredNoise', 'trialOutcome', 'noiseTimesMS', 'noiseCohsPC', 'signalCohPC'};
  missing = requiredFields(~isfield(sessionNoise, requiredFields));
  if ~isempty(missing)
    error('makeBetaKernel:MissingFields', '%s is missing: %s', fileInfo{iSession, 'fileName'}, strjoin(missing, ', '));
  end

  frameRateHz = headerScalar(H, 'frameRateHz');
  preStepMS = headerScalar(H, 'preStepMS');
  stepMS = headerScalar(H, 'stepMS');
  prefCohNoisePC = headerScalar(H, 'prefCohNoisePC');

  msPerVFrame = 1000 / frameRateHz;
  nFrames = round((preStepMS + stepMS) / msPerVFrame);
  tMS = (0:nFrames-1) * msPerVFrame;

  if ~initialized
    firstFrameRateHz = frameRateHz;
    firstPreStepMS = preStepMS;
    firstStepMS = stepMS;
    firstNFrames = nFrames;

    sessionKernelsPC = nan(nSessions, nFrames);
    sessionKernelSEMPC = nan(nSessions, nFrames);
    sessionMeanCorrectPC = nan(nSessions, nFrames);
    sessionMeanErrorPC = nan(nSessions, nFrames);
    initialized = true;
  else
    if frameRateHz ~= firstFrameRateHz || ...
        preStepMS ~= firstPreStepMS || ...
        stepMS ~= firstStepMS || ...
        nFrames ~= firstNFrames
      error('makeBetaKernel:IncompatibleSessions', '%s has incompatible timing or frame rate.', ...
              fileInfo{iSession, 'fileName'});
    end
  end

  useTrials = logical(sessionNoise.hasPreferredNoise(:));
  sessionNAllIncrementTrials(iSession) = sessionNoise.nTrials;
  sessionNTrials(iSession) = sum(useTrials);

  if sessionNTrials(iSession) == 0
    error('makeBetaKernel:NoPreferredNoiseTrials', '%s contains no preferred-noise trials.', ...
              fileInfo{iSession, 'fileName'});
  end

  selectedOutcome = double(sessionNoise.trialOutcome(useTrials));
  selectedOutcome = selectedOutcome(:);
  selectedTimes = sessionNoise.noiseTimesMS(useTrials);
  selectedNoise = sessionNoise.noiseCohsPC(useTrials);
  selectedSignal = double(sessionNoise.signalCohPC(useTrials));

  noiseMatrix = nan(sessionNTrials(iSession), nFrames);

  for iTrial = 1:sessionNTrials(iSession)
    times = double(selectedTimes{iTrial}(:)');
    values = double(selectedNoise{iTrial}(:)');

    if numel(times) ~= numel(values)
      error('makeBetaKernel:NoiseLengthMismatch', ...
        ['Preferred-noise trial %d in %s has %d times and %d values. ' ...
        'Length mismatches are tolerated only for no-preferred-noise trials.'], ...
        iTrial, fileInfo{iSession, 'fileName'}, numel(times), numel(values));
    end
    noiseMatrix(iTrial, :) = samplePiecewiseConstant(times, values, tMS);
  end

  isCorrect = selectedOutcome == 0;
  isError = selectedOutcome == 1;

  sessionNCorrect(iSession) = sum(isCorrect);
  sessionNError(iSession) = sum(isError);
  sessionSignalCohPC{iSession} = unique(selectedSignal(:)');
  sessionPrefCohNoisePC(iSession) = prefCohNoisePC;

  if sessionNCorrect(iSession) == 0 || sessionNError(iSession) == 0
    error('makeBetaKernel:EmptyOutcomeClass', '%s has %d correct and %d error preferred-noise trials.', ...
      fileInfo{iSession, 'fileName'}, sessionNCorrect(iSession), sessionNError(iSession));
  end
  [sessionKernel, sessionSEM, meanCorrect, meanError] = correctMinusErrorKernel(noiseMatrix, isCorrect, isError);

  sessionKernelsPC(iSession, :) = sessionKernel;
  sessionKernelSEMPC(iSession, :) = sessionSEM; %#ok<*AGROW>
  sessionMeanCorrectPC(iSession, :) = meanCorrect;
  sessionMeanErrorPC(iSession, :) = meanError;

  allNoise = [allNoise; noiseMatrix];
  allOutcome = [allOutcome; selectedOutcome];
  allSessionIndex = [allSessionIndex; repmat(iSession, sessionNTrials(iSession), 1)];
end
isCorrect = allOutcome == 0;
isError = allOutcome == 1;

[kernelPC, kernelSEMPC, meanCorrectPC, meanErrorPC] = correctMinusErrorKernel(allNoise, isCorrect, isError);

kernelData = struct();
kernelData.version = 2;
kernelData.method = ['trial-pooled change-side preferred-noise mean(hit) minus mean(miss); ' ...
                     'uses only sessionNoise.hasPreferredNoise trials'];
kernelData.timeOrigin = 'preStep onset';
kernelData.stepOnsetMS = firstPreStepMS;
kernelData.tMS = tMS;
kernelData.frameRateHz = firstFrameRateHz;
kernelData.msPerVFrame = 1000 / firstFrameRateHz;
kernelData.preStepMS = firstPreStepMS;
kernelData.stepMS = firstStepMS;

kernelData.kernelPC = kernelPC;
kernelData.kernelSEMPC = kernelSEMPC;
kernelData.meanCorrectPC = meanCorrectPC;
kernelData.meanErrorPC = meanErrorPC;
kernelData.nCorrect = sum(isCorrect);
kernelData.nError = sum(isError);
kernelData.nTrials = numel(allOutcome);

kernelData.sessionFileNames = fileInfo.fileName;
kernelData.sessionKernelsPC = sessionKernelsPC;
kernelData.sessionKernelSEMPC = sessionKernelSEMPC;
kernelData.sessionMeanCorrectPC = sessionMeanCorrectPC;
kernelData.sessionMeanErrorPC = sessionMeanErrorPC;

kernelData.sessionNCorrect = sessionNCorrect;
kernelData.sessionNError = sessionNError;
kernelData.sessionNTrials = sessionNTrials;
kernelData.sessionNAllIncrementTrials = sessionNAllIncrementTrials;
kernelData.sessionSignalCohPC = sessionSignalCohPC;
kernelData.sessionPrefCohNoisePC = sessionPrefCohNoisePC;

kernelData.trialSessionIndex = allSessionIndex;
kernelData.createdBy = mfilename;
kernelData.createdDate = datetime('now');

outputPath = fullfile(outputFolder, sprintf('BetaKernel_%s.mat', R.Animal));
save(outputPath, 'kernelData', '-v7.3');

makeBetaWeights(outputFolder, R.Animal);

% Plot the kernel
fig = figure('WindowStyle', 'docked');
hold on;
ci95 = 1.96 * kernelData.kernelSEMPC;
fill([kernelData.tMS, fliplr(kernelData.tMS)], [kernelData.kernelPC + ci95, fliplr(kernelData.kernelPC - ci95)], ...
  [0.85 0.85 1.00], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(kernelData.tMS, kernelData.kernelPC, 'b-', 'LineWidth', 1.5);
xline(kernelData.stepOnsetMS, 'k--', 'Step onset');
yline(0, 'k:');
xlabel('Time from preStep onset (ms)');
ylabel('Kernel (% coherence)');
title(sprintf('Preferred Direction Kernel (%d sessions, %d trials)', nSessions, kernelData.nTrials));
box off
exportgraphics(fig, fullfile(plotFolder, sprintf('AverageKernel_Preferred_%s.pdf', R.Animal)), ...
    'ContentType', 'vector');

end

%%-------------------------------------
function weightData = makeBetaWeights(dataFolder, animal)
% makeBetaWeights  Create normalized all-session and leave-one-out weights.
%
% Reads:
%   Data/AcrossOffsetSummaries/BetaKernel_AllSessions.mat
%
% Saves:
%   Data/AcrossOffsetSummaries/BetaWeights.mat
%
% The weighting function is the unsmoothed preferred-noise kernel over the
% 250-ms step interval. Weights are normalized to sum to 1, so a constant
% +1% coherence-noise waveform produces +1% effective coherence.
%
% Leave-one-session-out kernels are reconstructed exactly from the saved
% per-session correct/error means and trial counts.

% baseFolder = domainFolder(mfilename('fullpath'));
% acrossFolder = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries');
kernelPath = fullfile(dataFolder, sprintf('BetaKernel_%s.mat', animal));
outputPath = fullfile(dataFolder, sprintf('BetaWeights_%s.mat', animal));
S = load(kernelPath, 'kernelData');
K = S.kernelData;
requiredFields = { 'tMS', 'stepOnsetMS', 'stepMS', 'sessionFileNames', 'sessionMeanCorrectPC', ...
  'sessionMeanErrorPC', 'sessionNCorrect', 'sessionNError'};
missing = requiredFields(~isfield(K, requiredFields));
if ~isempty(missing)
  error('makeBetaWeights:MissingFields', ...
    'kernelData is missing: %s', strjoin(missing, ', '));
end
tMS = double(K.tMS(:)');
stepStartMS = double(K.stepOnsetMS);
stepEndMS = stepStartMS + double(K.stepMS);
stepMask = tMS >= stepStartMS & tMS < stepEndMS;
stepTMS = tMS(stepMask);
nSessions = numel(K.sessionFileNames);
meanCorrect = K.sessionMeanCorrectPC;
meanError = K.sessionMeanErrorPC;
nCorrect = double(K.sessionNCorrect(:));
nError = double(K.sessionNError(:));

if size(meanCorrect,1) ~= nSessions || size(meanError,1) ~= nSessions
  error('makeBetaWeights:SessionCountMismatch', 'Per-session mean arrays do not match sessionFileNames.');
end
if numel(nCorrect) ~= nSessions || numel(nError) ~= nSessions
  error('makeBetaWeights:CountMismatch', 'Per-session trial counts do not match sessionFileNames.');
end
if any(nCorrect <= 0) || any(nError <= 0)
  error('makeBetaWeights:EmptyOutcomeClass', 'Every session must contain both correct and error trials.');
end

% ---- All-session pooled kernel ----
allKernelPC = pooledKernel(meanCorrect, meanError, nCorrect, nError, true(nSessions,1));
allStepKernelPC = allKernelPC(stepMask);
allKernelSum = sum(allStepKernelPC);

assertUsableKernelSum(allKernelSum, 'all-session');

allWeights = allStepKernelPC / allKernelSum;

% ---- Leave-one-session-out kernels and weights ----
looKernelPC = nan(nSessions, numel(tMS));
looStepKernelPC = nan(nSessions, numel(stepTMS));
looWeights = nan(nSessions, numel(stepTMS));
looKernelSum = nan(nSessions, 1);
looNCorrect = nan(nSessions, 1);
looNError = nan(nSessions, 1);

for iSession = 1:nSessions
  includeMask = true(nSessions,1);
  includeMask(iSession) = false;

  thisKernel = pooledKernel(meanCorrect, meanError, ...
    nCorrect, nError, includeMask);
  thisStepKernel = thisKernel(stepMask);
  thisSum = sum(thisStepKernel);

  assertUsableKernelSum(thisSum, sprintf('leave-one-out session %d (%s)', iSession, K.sessionFileNames{iSession}));

  looKernelPC(iSession,:) = thisKernel;
  looStepKernelPC(iSession,:) = thisStepKernel;
  looWeights(iSession,:) = thisStepKernel / thisSum;
  looKernelSum(iSession) = thisSum;
  looNCorrect(iSession) = sum(nCorrect(includeMask));
  looNError(iSession) = sum(nError(includeMask));
end

% ---- Output structure ----
weightData = struct();
weightData.version = 1;
weightData.method = ...
  ['unsmoothed preferred-noise kernel over the step interval, ' ...
   'normalized so weights sum to 1'];

weightData.timeOrigin = K.timeOrigin;
weightData.stepOnsetMS = stepStartMS;
weightData.stepMS = K.stepMS;
weightData.tMS = tMS;
weightData.stepMask = stepMask;
weightData.stepTMS = stepTMS;

weightData.allSessionKernelPC = allKernelPC;
weightData.allSessionStepKernelPC = allStepKernelPC;
weightData.allSessionKernelSum = allKernelSum;
weightData.allSessionWeights = allWeights;
weightData.allSessionWeightSum = sum(allWeights);

weightData.sessionFileNames = K.sessionFileNames;
weightData.leaveOneOutKernelPC = looKernelPC;
weightData.leaveOneOutStepKernelPC = looStepKernelPC;
weightData.leaveOneOutKernelSum = looKernelSum;
weightData.leaveOneOutWeights = looWeights;
weightData.leaveOneOutWeightSums = sum(looWeights, 2);
weightData.leaveOneOutNCorrect = looNCorrect;
weightData.leaveOneOutNError = looNError;

weightData.nSessions = nSessions;
weightData.createdBy = mfilename;
weightData.createdDate = datetime('now');

% ---- Final consistency checks ----
if abs(weightData.allSessionWeightSum - 1) > 1e-12
  error('makeBetaWeights:AllWeightNormalizationFailed', 'All-session weights do not sum to 1.');
end

if any(abs(weightData.leaveOneOutWeightSums - 1) > 1e-12)
  error('makeBetaWeights:LOOWeightNormalizationFailed', 'At least one leave-one-out weight vector does not sum to 1.');
end

save(outputPath, 'weightData', '-v7.3');
end

% -------------------------------------------------------------------------
function kernelPC = pooledKernel(meanCorrect, meanError, nCorrect, nError, includeMask)

totalCorrect = sum(nCorrect(includeMask));
totalError = sum(nError(includeMask));
if totalCorrect <= 0 || totalError <= 0
  error('makeBetaWeights:NoOutcomeTrials', 'Pooled data must contain both correct and error trials.');
end
pooledCorrect = sum(meanCorrect(includeMask,:) .* nCorrect(includeMask), 1) / totalCorrect;
pooledError = sum(meanError(includeMask,:) .* nError(includeMask), 1) / totalError;
kernelPC = pooledCorrect - pooledError;
end

% -------------------------------------------------------------------------
function assertUsableKernelSum(kernelSum, label)

if ~isfinite(kernelSum)
  error('makeBetaWeights:NonfiniteKernelSum', 'The %s kernel sum is nonfinite.', label);
end

if kernelSum <= 0
  error('makeBetaWeights:NonpositiveKernelSum', 'The %s kernel sum is %.6g; expected a positive value.', ...
    label, kernelSum);
end

% Guard against pathological normalization.
if abs(kernelSum) < 1e-9
  error('makeBetaWeights:NearZeroKernelSum', 'The %s kernel sum is too close to zero for stable normalization.', label);
end
end


% -------------------------------------------------------------------------
function sampled = samplePiecewiseConstant(timesMS, valuesPC, queryTimesMS)

if isempty(timesMS) || isempty(valuesPC)
  error('makeBetaKernel:EmptyNoiseSequence', 'Preferred-noise trials must contain nonempty times and values.');
end
if numel(timesMS) ~= numel(valuesPC)
  error('makeBetaKernel:NoiseLengthMismatch','Noise times and values must have equal length.');
end
if any(~isfinite(timesMS)) || any(~isfinite(valuesPC))
  error('makeBetaKernel:NonfiniteNoise', 'Noise times and values must be finite.');
end
if any(diff(timesMS) < 0)
  error('makeBetaKernel:DecreasingTimes', 'Noise times must not decrease.');
end
if timesMS(1) ~= 0
  error('makeBetaKernel:BadStartTime', 'Noise times must begin at t = 0.');
end
sampled = nan(size(queryTimesMS));
eventIndex = 1;
currentValue = valuesPC(1);

for iQuery = 1:numel(queryTimesMS)
  queryTime = queryTimesMS(iQuery);

  while eventIndex < numel(timesMS) && ...
      timesMS(eventIndex + 1) <= queryTime
    eventIndex = eventIndex + 1;
    currentValue = valuesPC(eventIndex);
  end

  sampled(iQuery) = currentValue;
end
end

% -------------------------------------------------------------------------
function [kernel, kernelSEM, meanCorrect, meanError] = correctMinusErrorKernel(noiseMatrix, isCorrect, isError)

meanCorrect = mean(noiseMatrix(isCorrect, :), 1);
meanError = mean(noiseMatrix(isError, :), 1);
kernel = meanCorrect - meanError;

varCorrect = var(noiseMatrix(isCorrect, :), 0, 1);
varError = var(noiseMatrix(isError, :), 0, 1);

nCorrect = sum(isCorrect);
nError = sum(isError);

kernelSEM = sqrt(varCorrect / nCorrect + varError / nError);
end

% -------------------------------------------------------------------------
function value = headerScalar(H, fieldName)

if ~isfield(H, fieldName)
  error('makeBetaKernel:MissingHeaderField', 'sessionHeader.%s is required.', fieldName);
end

value = H.(fieldName);

while isstruct(value) && isfield(value, 'data')
  value = value.data;
end

if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  error('makeBetaKernel:BadHeaderField', 'sessionHeader.%s must be a numeric scalar.', fieldName);
end

value = double(value);
end