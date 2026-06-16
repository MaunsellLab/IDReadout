function kernelData = makeBetaKernel(omitSession)
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

if nargin < 1
  omitSession = [];
end

baseFolder = domainFolder(mfilename('fullpath'));
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
outputFolder = validFolder(fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions'));

files = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(files)
  error('makeBetaKernel:NoSessionData', ...
    'No session data files found in %s.', sessionDataFolder);
end

[~, order] = sort({files.name});
files = files(order);
nSessions = numel(files);

omitIndex = resolveOmittedSession(omitSession, files, sessionDataFolder);
includedSessionMask = true(1, nSessions);
if ~isempty(omitIndex)
  includedSessionMask(omitIndex) = false;
end

if ~any(includedSessionMask)
  error('makeBetaKernel:NoIncludedSessions',  'The requested omission leaves no sessions to analyze.');
end

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
  filePath = fullfile(sessionDataFolder, files(iSession).name);
  S = load(filePath, 'sessionNoise');

  if ~isfield(S, 'sessionNoise')
    error('makeBetaKernel:MissingSessionNoise', '%s does not contain sessionNoise.', files(iSession).name);
  end

  N = S.sessionNoise;
  H = N.sessionHeader;

  requiredFields = {'hasPreferredNoise', 'trialOutcome', 'noiseTimesMS', 'noiseCohsPC', 'signalCohPC'};
  missing = requiredFields(~isfield(N, requiredFields));
  if ~isempty(missing)
    error('makeBetaKernel:MissingFields', ...
      '%s is missing: %s', files(iSession).name, strjoin(missing, ', '));
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
        '%s has incompatible timing or frame rate.', files(iSession).name);
    end
  end

  useTrials = logical(N.hasPreferredNoise(:));
  sessionNAllIncrementTrials(iSession) = N.nTrials;
  sessionNTrials(iSession) = sum(useTrials);

  if sessionNTrials(iSession) == 0
    error('makeBetaKernel:NoPreferredNoiseTrials', ...
      '%s contains no preferred-noise trials.', files(iSession).name);
  end

selectedOutcome = double(N.trialOutcome(useTrials));
selectedOutcome = selectedOutcome(:);
selectedTimes = N.noiseTimesMS(useTrials);
  selectedNoise = N.noiseCohsPC(useTrials);
  selectedSignal = double(N.signalCohPC(useTrials));

  noiseMatrix = nan(sessionNTrials(iSession), nFrames);

  for iTrial = 1:sessionNTrials(iSession)
    times = double(selectedTimes{iTrial}(:)');
    values = double(selectedNoise{iTrial}(:)');

    if numel(times) ~= numel(values)
      error('makeBetaKernel:NoiseLengthMismatch', ...
        ['Preferred-noise trial %d in %s has %d times and %d values. ' ...
         'Length mismatches are tolerated only for no-preferred-noise trials.'], ...
        iTrial, files(iSession).name, numel(times), numel(values));
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
    error('makeBetaKernel:EmptyOutcomeClass', ...
      '%s has %d correct and %d error preferred-noise trials.', ...
      files(iSession).name, sessionNCorrect(iSession), sessionNError(iSession));
  end

  [sessionKernel, sessionSEM, meanCorrect, meanError] = ...
    correctMinusErrorKernel(noiseMatrix, isCorrect, isError);

  sessionKernelsPC(iSession, :) = sessionKernel;
  sessionKernelSEMPC(iSession, :) = sessionSEM;
  sessionMeanCorrectPC(iSession, :) = meanCorrect;
  sessionMeanErrorPC(iSession, :) = meanError;

  if includedSessionMask(iSession)
    allNoise = [allNoise; noiseMatrix]; %#ok<AGROW>
    allOutcome = [allOutcome; selectedOutcome]; %#ok<AGROW>
    allSessionIndex = [allSessionIndex; ...
      repmat(iSession, sessionNTrials(iSession), 1)]; %#ok<AGROW>
  end
end

isCorrect = allOutcome == 0;
isError = allOutcome == 1;

[kernelPC, kernelSEMPC, meanCorrectPC, meanErrorPC] = ...
  correctMinusErrorKernel(allNoise, isCorrect, isError);

kernelData = struct();
kernelData.version = 2;
kernelData.method = ...
  ['trial-pooled change-side preferred-noise mean(correct) minus mean(error); ' ...
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

kernelData.sessionFileNames = {files.name};
kernelData.includedSessionMask = includedSessionMask;
kernelData.omittedSession = omittedSessionDescription(omitIndex, files);

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

if isempty(omitIndex)
  outputName = 'BetaKernel_AllSessions.mat';
else
  [~, baseName] = fileparts(files(omitIndex).name);
  outputName = sprintf('BetaKernel_LeaveOut_%s.mat', baseName);
end

outputPath = fullfile(outputFolder, outputName);
save(outputPath, 'kernelData', '-v7.3');

% fprintf('Saved %s\n', outputPath);
% fprintf('Included %d sessions and %d preferred-noise trials: %d correct, %d error.\n', ...
%   sum(includedSessionMask), kernelData.nTrials, kernelData.nCorrect, kernelData.nError);

plotBetaKernel(kernelData);
end

% -------------------------------------------------------------------------
function sampled = samplePiecewiseConstant(timesMS, valuesPC, queryTimesMS)

if isempty(timesMS) || isempty(valuesPC)
  error('makeBetaKernel:EmptyNoiseSequence', ...
    'Preferred-noise trials must contain nonempty times and values.');
end

if numel(timesMS) ~= numel(valuesPC)
  error('makeBetaKernel:NoiseLengthMismatch', ...
    'Noise times and values must have equal length.');
end

if any(~isfinite(timesMS)) || any(~isfinite(valuesPC))
  error('makeBetaKernel:NonfiniteNoise', ...
    'Noise times and values must be finite.');
end

if any(diff(timesMS) < 0)
  error('makeBetaKernel:DecreasingTimes', ...
    'Noise times must not decrease.');
end

if timesMS(1) ~= 0
  error('makeBetaKernel:BadStartTime', ...
    'Noise times must begin at t = 0.');
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
function [kernel, kernelSEM, meanCorrect, meanError] = ...
  correctMinusErrorKernel(noiseMatrix, isCorrect, isError)

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
function omitIndex = resolveOmittedSession(omitSession, files, dataFolder)

omitIndex = [];

if isempty(omitSession)
  return;
end

if isnumeric(omitSession)
  validateattributes(omitSession, {'numeric'}, ...
    {'scalar', 'integer', 'positive', '<=', numel(files)});
  omitIndex = double(omitSession);
  return;
end

if ~(ischar(omitSession) || isstring(omitSession))
  error('makeBetaKernel:BadOmitSession', ...
    'omitSession must be empty, a numeric index, or a file name.');
end

query = char(omitSession);
queryBase = erase(query, '.mat');
matches = false(1, numel(files));

for iFile = 1:numel(files)
  [~, fileBase] = fileparts(files(iFile).name);
  matches(iFile) = strcmpi(query, files(iFile).name) || ...
    strcmpi(queryBase, fileBase);

  if ~matches(iFile)
    S = load(fullfile(dataFolder, files(iFile).name), 'sessionNoise');
    if isfield(S, 'sessionNoise') && isfield(S.sessionNoise, 'sourceFile')
      [~, sourceBase] = fileparts(S.sessionNoise.sourceFile);
      matches(iFile) = strcmpi(query, S.sessionNoise.sourceFile) || ...
        strcmpi(queryBase, sourceBase);
    end
  end
end

if sum(matches) ~= 1
  error('makeBetaKernel:OmitSessionNotUnique', ...
    'omitSession matched %d files; expected exactly one.', sum(matches));
end

omitIndex = find(matches);
end

% -------------------------------------------------------------------------
function description = omittedSessionDescription(omitIndex, files)

if isempty(omitIndex)
  description = '';
else
  description = files(omitIndex).name;
end
end

% -------------------------------------------------------------------------
function value = headerScalar(H, fieldName)

if ~isfield(H, fieldName)
  error('makeBetaKernel:MissingHeaderField', ...
    'sessionHeader.%s is required.', fieldName);
end

value = H.(fieldName);

while isstruct(value) && isfield(value, 'data')
  value = value.data;
end

if ~(isnumeric(value) || islogical(value)) || ~isscalar(value)
  error('makeBetaKernel:BadHeaderField', ...
    'sessionHeader.%s must be a numeric scalar.', fieldName);
end

value = double(value);
end

% -------------------------------------------------------------------------
function plotBetaKernel(kernelData)

figure;
hold on

ci95 = 1.96 * kernelData.kernelSEMPC;

fill([kernelData.tMS, fliplr(kernelData.tMS)], ...
     [kernelData.kernelPC + ci95, ...
      fliplr(kernelData.kernelPC - ci95)], ...
     [0.85 0.85 1.00], ...
     'EdgeColor', 'none', ...
     'FaceAlpha', 0.4);

plot(kernelData.tMS, kernelData.kernelPC, ...
  'b-', 'LineWidth', 1.5);

xline(kernelData.stepOnsetMS, 'k--', 'Step onset');
yline(0, 'k:');

xlabel('Time from preStep onset (ms)');
ylabel('Kernel (% coherence)');
title(sprintf('Preferred-noise kernel: %d sessions, %d trials', ...
  sum(kernelData.includedSessionMask), kernelData.nTrials));

box off
end
