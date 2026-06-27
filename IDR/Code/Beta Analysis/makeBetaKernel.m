function kernelData = makeBetaKernel()
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

baseFolder = domainFolder(mfilename('fullpath'));
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', 'BetaAnalysis');
outputFolder = validFolder(fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries'));
plotFolder = validFolder(fullfile(baseFolder, 'Plots', 'AcrossProbes', 'Kernels'));


[selectedFiles, fileInfo] = selectAnalysisFiles(sessionDataFolder);
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
    error('makeBetaKernel:MissingFields', ...
      '%s is missing: %s', fileInfo{iSession, 'fileName'}, strjoin(missing, ', '));
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
      error('makeBetaKernel:IncompatibleSessions', ...
        '%s has incompatible timing or frame rate.', fileInfo{iSession, 'fileName'});
    end
  end

  useTrials = logical(sessionNoise.hasPreferredNoise(:));
  sessionNAllIncrementTrials(iSession) = sessionNoise.nTrials;
  sessionNTrials(iSession) = sum(useTrials);

  if sessionNTrials(iSession) == 0
    error('makeBetaKernel:NoPreferredNoiseTrials', ...
      '%s contains no preferred-noise trials.', fileInfo{iSession, 'fileName'});
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

    noiseMatrix(iTrial, :) = ...
      samplePiecewiseConstant(times, values, tMS);
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

  [sessionKernel, sessionSEM, meanCorrect, meanError] = ...
    correctMinusErrorKernel(noiseMatrix, isCorrect, isError);

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
kernelData.method = ['trial-pooled change-side preferred-noise mean(correct) minus mean(error); ' ...
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

outputName = 'BetaKernel_AllSessions.mat';
outputPath = fullfile(outputFolder, outputName);
save(outputPath, 'kernelData', '-v7.3');

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
exportgraphics(fig, fullfile(plotFolder, 'AverageKernel_Preferred.pdf'), 'ContentType','vector');

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