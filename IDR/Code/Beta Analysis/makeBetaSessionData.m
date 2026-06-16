function [nCreated, nSkipped] = makeBetaSessionData(replace)
% makeBetaSessionData  Extract full-session increment trials for beta analyses.
%
% Creates one file per full recording session in:
%   Data/FullSessions/BetaAnalysis/SessionData
%
% Each file contains sessionNoise with all valid coherence-increment trials,
% including both noise and no-noise trials and all tested signal coherences.
% Preferred-direction noise is taken from the change-side patch and retained
% in its original piecewise-constant representation.
%
% Usage:
%   makeBetaSessionData();       % create missing files
%   makeBetaSessionData(true);   % replace existing files
%
% OUTPUT
%   nCreated  Number of session files written
%   nSkipped  Number of existing session files skipped

% cleanupObj = initProjectPath(); %#ok<NASGU>

if nargin < 1 || isempty(replace)
  replace = false;
end
validateattributes(replace, {'logical'}, {'scalar'});

% sessionFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions');
% outputFolder = validFolder(fullfile(sessionFolder, ...
%   'BetaAnalysis', 'SessionData'));
% 
% files = dir(fullfile(sessionFolder, '*.mat'));
% files = files(~endsWith({files.name}, '_fileInfo.mat'));
% 
% if isempty(files)
%   fprintf('No full-session .mat files found in %s\n', sessionFolder);
%   nCreated = 0;
%   nSkipped = 0;
%   return;
% end

sessionFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions');
betaFolder = fullfile(sessionFolder, 'BetaAnalysis');
outputFolder = fullfile(betaFolder, 'SessionData');
acrossFolder = fullfile(betaFolder, 'AcrossSessions');

if replace
  if isfolder(outputFolder)
    rmdir(outputFolder, 's');
  end
  if isfolder(acrossFolder)
    rmdir(acrossFolder, 's');
  end
end

validFolder(outputFolder);
validFolder(acrossFolder);
[selectedFiles, fileInfo] = selectAnalysisFiles(sessionFolder);
if isempty(selectedFiles)
  error('makeBetaSessionData:NoSelectedFiles', 'No full-session files passed selectAnalysisFiles.');
end

nCreated = 0;
nSkipped = 0;

for iFile = 1:numel(selectedFiles)

    inputPath = selectedFiles{iFile};
    [~, baseName] = fileparts(inputPath);
    outputPath = fullfile(outputFolder, [baseName '.mat']);

    if isfile(outputPath) && ~replace
        nSkipped = nSkipped + 1;
        continue
    end
    nCreated = nCreated + 1;

    fprintf('      processing %s ...\n', baseName);

    S = load(inputPath, 'sessionHeader', 'trials');

    sessionNoise = extractSessionNoise( ...
        S.sessionHeader, S.trials);

    [~, sourceName, sourceExt] = fileparts(inputPath);
    sessionNoise.sourceFile = [sourceName sourceExt];

    save(outputPath, 'sessionNoise', '-v7.3');
end

selectionManifest = struct();
selectionManifest.fileInfo = fileInfo;
selectionManifest.selectedFiles = selectedFiles;
selectionManifest.createdBy = mfilename;
selectionManifest.createdDate = datetime('now');

save(fullfile(acrossFolder, 'BetaSelectionManifest.mat'), 'selectionManifest');
% fprintf('Created %d files; skipped %d existing files.\n', ...
%   nCreated, nSkipped);
end

% -------------------------------------------------------------------------
function sessionNoise = extractSessionNoise(sessionHeader, trials)
% Extract all certified increment trials with valid behavioral outcomes.

requiredHeaderFields = {'preStepMS', 'stepMS', 'cohNoiseFrameMS'};
for iField = 1:numel(requiredHeaderFields)
  f = requiredHeaderFields{iField};
  assert(isfield(sessionHeader, f), ...
    'makeBetaSessionData:MissingSessionHeaderField', ...
    'sessionHeader.%s is required.', f);
end

% Acquisition constants expected throughout this analysis.
assertScalarEqual(sessionHeader.preStepMS, 750, 'preStepMS');
assertScalarEqual(sessionHeader.stepMS, 250, 'stepMS');
assertScalarEqual(sessionHeader.cohNoiseFrameMS, 50, ...
  'cohNoiseFrameMS');

nAllTrials = numel(trials);
useTrial = false(1, nAllTrials);

for iTrial = 1:nAllTrials
  tr = trials{iTrial};

  if ~hasRequiredTrialFields(tr)
    continue;
  end

  D = tr.trial.data;
  trialOutcome = scalarValue(tr.trialEnd);
  trialCertify = scalarValue(tr.trialCertify);

  % Raw changeIndex convention: 0=DEC, 1=INC.
  isIncrement = scalarValue(D.changeIndex) == 1;
  isValidOutcome = ismember(trialOutcome, [0 1]);
  isCertified = trialCertify == 0;

  useTrial(iTrial) = isIncrement && isValidOutcome && isCertified;
end

trialIdx = find(useTrial);
assert(~isempty(trialIdx), ...
  'makeBetaSessionData:NoTrials', ...
  'No valid coherence-increment trials were found.');

nTrials = numel(trialIdx);

trialOutcome = nan(1, nTrials);
changeSide = nan(1, nTrials);

baseCohPC = nan(1, nTrials);
stepCohPC = nan(1, nTrials);
signalCohPC = nan(1, nTrials);

hasCohNoise = false(1, nTrials);
prefCohNoisePC = nan(1, nTrials);
hasPreferredNoise = false(1, nTrials);

noiseTimesMS = cell(1, nTrials);
noiseCohsPC = cell(1, nTrials);

for k = 1:nTrials
  iTrial = trialIdx(k);
  tr = trials{iTrial};
  D = tr.trial.data;

  trialOutcome(k) = scalarValue(tr.trialEnd);
  changeSide(k) = scalarValue(D.changeSide);

  assert(ismember(changeSide(k), [0 1]), ...
    'makeBetaSessionData:BadChangeSide', ...
    'Trial %d has changeSide=%g; expected 0 (RF) or 1 (Opp).', ...
    iTrial, changeSide(k));

  baseCohPC(k) = scalarValue(D.baseCohPC);
  stepCohPC(k) = scalarValue(D.stepCohPC);
  signalCohPC(k) = stepCohPC(k) - baseCohPC(k);

  assert(signalCohPC(k) > 0, ...
    'makeBetaSessionData:NonIncrementCoherence', ...
    ['Trial %d is marked as an increment but has ' ...
     'stepCohPC=%g and baseCohPC=%g.'], ...
    iTrial, stepCohPC(k), baseCohPC(k));

  hasCohNoise(k) = logical(scalarValue(D.cohNoise));
  prefCohNoisePC(k) = scalarValue(D.prefCohNoisePC);

  % This is the flag relevant to the preferred-direction kernel and
  % effective-noise analyses.
  hasPreferredNoise(k) = ...
    hasCohNoise(k) && prefCohNoisePC(k) ~= 0;

  % change* variables describe the patch whose coherence stepped,
  % independent of whether that patch was RF or Opp.
  times = double(dataValue(tr.changeTimesMS));
  noise = double(dataValue(tr.changePrefCohsPC));

  times = times(:)';
  noise = noise(:)';

  assert(numel(times) <= numel(noise), ...
    'makeBetaSessionData:NoiseLengthMismatch', ...
    'Trial %d has %d noise times but %d preferred-noise values.', ...
    iTrial, numel(times), numel(noise));

  if isempty(times)
    % Empty arrays are acceptable only when no preferred noise was present.
    assert(~hasPreferredNoise(k), ...
      'makeBetaSessionData:MissingPreferredNoise', ...
      'Trial %d has preferred noise but no recorded noise sequence.', ...
      iTrial);
  else
    assert(times(1) == 0, ...
      'makeBetaSessionData:BadNoiseStart', ...
      'Trial %d noise times do not begin at 0 ms.', iTrial);

    dt = diff(times);
    assert(all(dt >= 0), ...
      'makeBetaSessionData:DecreasingNoiseTimes', ...
      'Trial %d noise times decrease.', iTrial);

    assert(all(times(dt == 0) == 0), ...
      'makeBetaSessionData:RepeatedNoiseTime', ...
      'Trial %d has a repeated noise time other than t = 0.', ...
      iTrial);

    assert(all(isfinite(times)) && all(isfinite(noise)), ...
      'makeBetaSessionData:NonfiniteNoiseSequence', ...
      'Trial %d has nonfinite noise times or values.', iTrial);
  end

  % A trial without preferred-direction noise should have no preferred
  % perturbation, even if the shared event-time array is populated because
  % another stream was active.
  if ~hasPreferredNoise(k) && ~isempty(noise)
    assert(all(noise == 0), ...
      'makeBetaSessionData:UnexpectedPreferredNoise', ...
      ['Trial %d is classified as having no preferred-direction noise, ' ...
       'but changePrefCohsPC contains nonzero values.'], iTrial);
  end

  noiseTimesMS{k} = times;
  noiseCohsPC{k} = noise;
end

sessionNoise = struct();
sessionNoise.version = 2;
sessionNoise.sessionHeader = sessionHeader;

sessionNoise.nTrials = nTrials;
sessionNoise.trialIdx = trialIdx;
sessionNoise.trialOutcome = trialOutcome;
sessionNoise.changeSide = changeSide;

sessionNoise.baseCohPC = baseCohPC;
sessionNoise.stepCohPC = stepCohPC;
sessionNoise.signalCohPC = signalCohPC;

sessionNoise.hasCohNoise = hasCohNoise;
sessionNoise.prefCohNoisePC = prefCohNoisePC;
sessionNoise.hasPreferredNoise = hasPreferredNoise;

sessionNoise.noiseTimesMS = noiseTimesMS;
sessionNoise.noiseCohsPC = noiseCohsPC;

sessionNoise.selection = struct( ...
  'stepType', 'increment', ...
  'rawChangeIndex', 1, ...
  'requiresCohNoise', false, ...
  'requiresNonzeroPrefCohNoisePC', false, ...
  'requiresTrialCertify', 0, ...
  'validTrialOutcomes', [0 1], ...
  'noisePatch', 'change side');

sessionNoise.counts = struct( ...
  'nTrials', nTrials, ...
  'nPreferredNoiseTrials', sum(hasPreferredNoise), ...
  'nNoPreferredNoiseTrials', sum(~hasPreferredNoise), ...
  'nSignalCoherences', numel(unique(signalCohPC)), ...
  'signalCoherencesPC', unique(signalCohPC));
end

% -------------------------------------------------------------------------
function tf = hasRequiredTrialFields(tr)

tf = isfield(tr, 'trial') && isfield(tr.trial, 'data') && ...
  isfield(tr, 'trialEnd') && isfield(tr, 'trialCertify') && ...
  isfield(tr, 'changeTimesMS') && ...
  isfield(tr, 'changePrefCohsPC');

if ~tf
  return;
end

D = tr.trial.data;
required = {'changeIndex', 'changeSide', 'cohNoise', ...
  'prefCohNoisePC', 'baseCohPC', 'stepCohPC'};

tf = all(isfield(D, required));
end

% -------------------------------------------------------------------------
function assertScalarEqual(value, expected, fieldName)

value = scalarValue(value);

assert(isfinite(value) && value == expected, ...
  'makeBetaSessionData:UnexpectedAcquisitionValue', ...
  'sessionHeader.%s is %g; expected %g.', ...
  fieldName, value, expected);
end

% -------------------------------------------------------------------------
function value = scalarValue(x)

value = dataValue(x);

assert(isnumeric(value) || islogical(value), ...
  'makeBetaSessionData:NonNumericScalar', ...
  'Expected a numeric or logical scalar.');

assert(isscalar(value), ...
  'makeBetaSessionData:NonScalarValue', ...
  'Expected a scalar value.');

value = double(value);
end

% -------------------------------------------------------------------------
function value = dataValue(x)

value = x;

while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
end
