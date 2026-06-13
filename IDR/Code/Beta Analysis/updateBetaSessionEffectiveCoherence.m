function updateBetaSessionEffectiveCoherence(replace)
% updateBetaSessionEffectiveCoherence
% Add trialwise effective-noise and effective-coherence values to each
% expanded BetaAnalysis SessionData file.
%
% Preferred-noise trials receive leave-one-out weighted noise.
% No-preferred-noise trials receive effectiveNoisePC = 0.
% All trials receive effectiveCohPC = signalCohPC + effectiveNoisePC.

cleanupObj = initProjectPath(); %#ok<NASGU>

if nargin < 1 || isempty(replace)
  replace = false;
end

baseFolder = folderPath();
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
weightPath = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions', 'BetaWeights.mat');

if ~isfile(weightPath)
  error('updateBetaSessionEffectiveCoherence:MissingWeights', ...
    'Weight file not found: %s', weightPath);
end

W = load(weightPath, 'weightData');
weightData = W.weightData;

files = dir(fullfile(sessionDataFolder, '*.mat'));
[~, order] = sort({files.name});
files = files(order);

weightNames = string(weightData.sessionFileNames(:));
stepTMS = double(weightData.stepTMS(:)');
nStepFrames = numel(stepTMS);

nUpdated = 0;
nSkipped = 0;

for iFile = 1:numel(files)
  fileName = files(iFile).name;
  filePath = fullfile(sessionDataFolder, fileName);

  match = find(strcmpi(weightNames, fileName));
  if numel(match) ~= 1
    error('updateBetaSessionEffectiveCoherence:SessionMatchFailure', ...
      '%s matched %d rows in weightData.sessionFileNames.', ...
      fileName, numel(match));
  end

  S = load(filePath, 'sessionNoise');
  sessionNoise = S.sessionNoise;

  if isfield(sessionNoise, 'trialAnalysis') && ~replace
    fprintf('Skipping %s: trialAnalysis already exists.\n', fileName);
    nSkipped = nSkipped + 1;
    continue;
  end

  nTrials = sessionNoise.nTrials;
  signalCohPC = double(sessionNoise.signalCohPC(:));
  hasPreferredNoise = logical(sessionNoise.hasPreferredNoise(:));

  if numel(signalCohPC) ~= nTrials || numel(hasPreferredNoise) ~= nTrials || ...
      numel(sessionNoise.noiseTimesMS) ~= nTrials || ...
      numel(sessionNoise.noiseCohsPC) ~= nTrials
    error('updateBetaSessionEffectiveCoherence:TrialCountMismatch', ...
      '%s has inconsistent trial-level array lengths.', fileName);
  end

  weights = double(weightData.leaveOneOutWeights(match, :));
  if numel(weights) ~= nStepFrames || abs(sum(weights) - 1) > 1e-12
    error('updateBetaSessionEffectiveCoherence:BadWeights', ...
      'Invalid leave-one-out weights for %s.', fileName);
  end

  stepNoisePC = zeros(nTrials, nStepFrames);
  effectiveNoisePC = zeros(nTrials, 1);
  meanStepNoisePC = zeros(nTrials, 1);

  preferredIdx = find(hasPreferredNoise);
  for j = 1:numel(preferredIdx)
    iTrial = preferredIdx(j);
    timesMS = double(sessionNoise.noiseTimesMS{iTrial}(:)');
    noisePC = double(sessionNoise.noiseCohsPC{iTrial}(:)');

    if numel(timesMS) ~= numel(noisePC)
      error('updateBetaSessionEffectiveCoherence:NoiseLengthMismatch', ...
        ['Preferred-noise trial %d in %s has %d times and %d values. ' ...
         'Length mismatches are tolerated only on no-preferred-noise trials.'], ...
        iTrial, fileName, numel(timesMS), numel(noisePC));
    end

    stepNoisePC(iTrial, :) = samplePiecewiseConstant(timesMS, noisePC, stepTMS);
    effectiveNoisePC(iTrial) = stepNoisePC(iTrial, :) * weights(:);
    meanStepNoisePC(iTrial) = mean(stepNoisePC(iTrial, :));
  end

  effectiveCohPC = signalCohPC + effectiveNoisePC;

  trialAnalysis = struct();
  trialAnalysis.version = 2;
  trialAnalysis.method = ...
    ['preferred-direction change-side noise sampled on step video frames; ' ...
     'weighted by leave-one-session-out kernel normalized to sum to 1; ' ...
     'no-preferred-noise trials assigned zero effective noise'];
  trialAnalysis.stepTMS = stepTMS;
  trialAnalysis.leaveOneOutWeights = weights;
  trialAnalysis.leaveOneOutKernelSum = weightData.leaveOneOutKernelSum(match);
  trialAnalysis.stepNoisePC = stepNoisePC;
  trialAnalysis.effectiveNoisePC = effectiveNoisePC;
  trialAnalysis.meanStepNoisePC = meanStepNoisePC;
  trialAnalysis.signalCohPC = signalCohPC;
  trialAnalysis.effectiveCohPC = effectiveCohPC;
  trialAnalysis.hasPreferredNoise = hasPreferredNoise;
  trialAnalysis.weightSourceFile = weightPath;
  trialAnalysis.omittedSessionFile = fileName;
  trialAnalysis.createdBy = mfilename;
  trialAnalysis.createdDate = datetime('now');

  sessionNoise.trialAnalysis = trialAnalysis;
  save(filePath, 'sessionNoise', '-append');

  fprintf(['Updated %s: %d total trials, %d preferred-noise, ' ...
           '%d no-preferred-noise.\n'], ...
    fileName, nTrials, sum(hasPreferredNoise), sum(~hasPreferredNoise));
  nUpdated = nUpdated + 1;
end

fprintf('Finished: %d files updated, %d skipped.\n', nUpdated, nSkipped);
end

function sampled = samplePiecewiseConstant(timesMS, valuesPC, queryTimesMS)
if isempty(timesMS) || isempty(valuesPC)
  error('updateBetaSessionEffectiveCoherence:EmptyNoiseSequence', ...
    'Preferred-noise trials must contain nonempty times and values.');
end
if numel(timesMS) ~= numel(valuesPC)
  error('updateBetaSessionEffectiveCoherence:NoiseLengthMismatch', ...
    'Noise times and values must have equal length.');
end
if any(~isfinite(timesMS)) || any(~isfinite(valuesPC))
  error('updateBetaSessionEffectiveCoherence:NonfiniteNoise', ...
    'Noise times and values must be finite.');
end
if any(diff(timesMS) < 0)
  error('updateBetaSessionEffectiveCoherence:DecreasingTimes', ...
    'Noise times must not decrease.');
end
if timesMS(1) ~= 0
  error('updateBetaSessionEffectiveCoherence:BadStartTime', ...
    'Noise times must begin at t = 0.');
end

sampled = nan(size(queryTimesMS));
eventIndex = 1;
currentValue = valuesPC(1);
for iQuery = 1:numel(queryTimesMS)
  queryTime = queryTimesMS(iQuery);
  while eventIndex < numel(timesMS) && timesMS(eventIndex + 1) <= queryTime
    eventIndex = eventIndex + 1;
    currentValue = valuesPC(eventIndex);
  end
  sampled(iQuery) = currentValue;
end
end
