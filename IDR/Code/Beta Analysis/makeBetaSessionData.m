function [nCreated, nSkipped] = makeBetaSessionData(varargin)
% makeBetaSessionData  Extract full-session increment trials for beta analyses.
%
% Creates one file per full recording session in:
%   Data/FullSessions/BetaAnalysis
%
% Each file contains sessionNoise with all valid coherence-increment trials,
% including both noise and no-noise trials and all tested signal coherences.
% Preferred-direction noise is taken from the change-side patch and retained
% in its original piecewise-constant representation.
%
% OUTPUT
%   nCreated  Number of session files written
%   nSkipped  Number of existing session files skipped

p = inputParser;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) isempty(x) || islogical(x));
parse(p, varargin{:});
opts = p.Results;

acrossFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries');
sessionFolder = char(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));
outputFolder = fullfile(sessionFolder, 'BetaAnalysis');

if opts.Replace
  if isfolder(outputFolder)
    rmdir(outputFolder, 's');
  end
  if isfolder(acrossFolder)
    rmdir(acrossFolder, 's');
  end
end

validFolder(outputFolder);
validFolder(acrossFolder);
[selectedFiles, ~] = selectAnalysisFiles(sessionFolder, 'Animal', opts.Animal);
if isempty(selectedFiles)
  error('makeBetaSessionData:NoSelectedFiles', 'No full-session files passed selectAnalysisFiles.');
end
nCreated = 0;
nSkipped = 0;
for iFile = 1:numel(selectedFiles)
  inputPath = selectedFiles{iFile};
  [~, baseName] = fileparts(inputPath);
  outputPath = fullfile(outputFolder, [baseName '.mat']);
  if isfile(outputPath) && ~opts.Replace
    nSkipped = nSkipped + 1;
    continue
  end
  nCreated = nCreated + 1;
  fprintf('      processing %s ...\n', baseName);
  load(inputPath, 'sessionHeader', 'trials');
  sessionNoise = extractSessionNoise(sessionHeader, trials);
  [~, sourceName, sourceExt] = fileparts(inputPath);
  sessionNoise.sourceFile = [sourceName sourceExt];
  save(outputPath, 'sessionHeader', 'sessionNoise', '-v7.3');
end
validateBetaSessionData();
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

%% -------------------------------------------------------------------------
function value = dataValue(x)

value = x;

while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
end

%% -------------------------------------------------------------------------
function [sessionTable, issueTable] = validateBetaSessionData()
% validateBetaSessionData  Validate expanded BetaAnalysis SessionData files.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/*.mat
%
% Supports sessionNoise.version >= 2, containing all valid increment trials,
% including both preferred-noise and no-preferred-noise trials.

dataFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions', 'BetaAnalysis');
files = dir(fullfile(dataFolder, '*.mat'));
if isempty(files)
  error('validateBetaSessionData:NoFiles', 'No SessionData files found in %s.', dataFolder);
end

[~, order] = sort({files.name});
files = files(order);

summaryRows = repmat(emptySummaryRow(), numel(files), 1);
issueRows = repmat(emptyIssueRow(), 0, 1);

for iFile = 1:numel(files)
  fileName = files(iFile).name;
  filePath = fullfile(dataFolder, fileName);

  row = emptySummaryRow();
  row.fileName = fileName;

  S = load(filePath, 'sessionNoise');
  if ~isfield(S, 'sessionNoise')
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      'File does not contain sessionNoise.'); %#ok<AGROW>
    row.nErrors = 1;
    row.nIssues = 1;
    summaryRows(iFile) = row;
    continue;
  end

  N = S.sessionNoise;
  requiredFields = {'version', 'sessionHeader', 'nTrials', 'trialIdx', 'trialOutcome', 'changeSide', 'baseCohPC', 
    'stepCohPC', 'signalCohPC', 'hasCohNoise', 'prefCohNoisePC', 'hasPreferredNoise', 'noiseTimesMS', 'noiseCohsPC'};

  missing = requiredFields(~isfield(N, requiredFields));
  if ~isempty(missing)
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      sprintf('Missing fields: %s', strjoin(missing, ', '))); %#ok<AGROW>
    row.nErrors = 1;
    row.nIssues = 1;
    summaryRows(iFile) = row;
    continue;
  end

  H = N.sessionHeader;
  row.version = double(N.version);
  row.nTrials = double(N.nTrials);
  row.preStepMS = getHeaderScalar(H, 'preStepMS');
  row.stepMS = getHeaderScalar(H, 'stepMS');
  row.cohNoiseFrameMS = getHeaderScalar(H, 'cohNoiseFrameMS');
  row.frameRateHz = getHeaderScalar(H, 'frameRateHz');
  row.headerPrefCohNoisePC = getHeaderScalar(H, 'prefCohNoisePC');

  n = N.nTrials;

  lengths = [numel(N.trialIdx), numel(N.trialOutcome), ...
    numel(N.changeSide), numel(N.baseCohPC), numel(N.stepCohPC), ...
    numel(N.signalCohPC), numel(N.hasCohNoise), ...
    numel(N.prefCohNoisePC), numel(N.hasPreferredNoise), ...
    numel(N.noiseTimesMS), numel(N.noiseCohsPC)];

  if any(lengths ~= n)
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      sprintf('nTrials=%d but trial-level lengths are %s.', ...
      n, mat2str(lengths))); %#ok<AGROW>
  end

  trialIdx = double(N.trialIdx(:));
  trialOutcome = double(N.trialOutcome(:));
  changeSide = double(N.changeSide(:));
  baseCohPC = double(N.baseCohPC(:));
  stepCohPC = double(N.stepCohPC(:));
  signalCohPC = double(N.signalCohPC(:));
  hasCohNoise = logical(N.hasCohNoise(:));
  prefCohNoisePC = double(N.prefCohNoisePC(:));
  hasPreferredNoise = logical(N.hasPreferredNoise(:));

  if any(~ismember(trialOutcome, [0 1]))
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      'trialOutcome contains values other than 0 or 1.'); %#ok<AGROW>
  end

  if any(~ismember(changeSide, [0 1]))
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      'changeSide contains values other than 0 or 1.'); %#ok<AGROW>
  end

  if any(trialIdx < 1) || any(mod(trialIdx,1) ~= 0) || ...
      numel(unique(trialIdx)) ~= numel(trialIdx)
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      'trialIdx must contain unique positive integers.'); %#ok<AGROW>
  end

  if any(~isfinite(baseCohPC)) || any(~isfinite(stepCohPC)) || ...
      any(~isfinite(signalCohPC))
    issueRows(end+1) = makeIssue(fileName, NaN, 'error', ...
      'Coherence vectors contain nonfinite values.'); %#ok<AGROW>
  end

  coherenceMismatch = abs(signalCohPC - (stepCohPC - baseCohPC)) > 1e-12;
  if any(coherenceMismatch)
    firstBad = find(coherenceMismatch, 1, 'first');
    issueRows(end+1) = makeIssue(fileName, trialIdx(firstBad), 'error', ...
      'signalCohPC does not equal stepCohPC - baseCohPC.'); %#ok<AGROW>
  end

  if any(signalCohPC <= 0)
    firstBad = find(signalCohPC <= 0, 1, 'first');
    issueRows(end+1) = makeIssue(fileName, trialIdx(firstBad), 'error', ...
      'Included increment trial has nonpositive signal coherence.'); %#ok<AGROW>
  end

  expectedHasPreferredNoise = hasCohNoise & prefCohNoisePC ~= 0;
  if any(hasPreferredNoise ~= expectedHasPreferredNoise)
    firstBad = find(hasPreferredNoise ~= expectedHasPreferredNoise, 1, 'first');
    issueRows(end+1) = makeIssue(fileName, trialIdx(firstBad), 'error', ...
      ['hasPreferredNoise does not equal hasCohNoise && ' ...
       'prefCohNoisePC ~= 0.']); %#ok<AGROW>
  end

  row.nPreferredNoiseTrials = sum(hasPreferredNoise);
  row.nNoPreferredNoiseTrials = sum(~hasPreferredNoise);
  row.nCohNoiseTrials = sum(hasCohNoise);
  row.nNoCohNoiseTrials = sum(~hasCohNoise);
  row.nSignalCoherences = numel(unique(signalCohPC));
  row.minSignalCohPC = min(signalCohPC);
  row.maxSignalCohPC = max(signalCohPC);
  row.fractionCorrect = mean(trialOutcome == 0);
  row.signalCoherencesPC = string(mat2str(unique(signalCohPC(:)')));

  allPositiveIntervals = [];
  finalIntervals = [];
  observedNoiseAbs = [];
  nDuplicateZeroTrials = 0;

  totalMS = row.preStepMS + row.stepMS;
  if isfinite(row.frameRateHz) && row.frameRateHz > 0
    timingToleranceMS = 1000 / row.frameRateHz / 2 + 1e-6;
  else
    timingToleranceMS = 1e-6;
  end

  for k = 1:n
    times = double(N.noiseTimesMS{k}(:)');
    values = double(N.noiseCohsPC{k}(:)');
    parentTrial = trialIdx(k);

    if isempty(times)
      if hasPreferredNoise(k)
        issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
          'Preferred-noise trial has an empty sequence.'); %#ok<AGROW>
      end
      continue;
    end

    if any(~isfinite(times)) || any(~isfinite(values))
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        'Noise times or values contain nonfinite entries.'); %#ok<AGROW>
      continue;
    end

    if times(1) ~= 0
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        sprintf('First recorded time is %.6g ms rather than 0.', times(1))); %#ok<AGROW>
    end

    dt = diff(times);

    if any(dt < 0)
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        'Noise times decrease.'); %#ok<AGROW>
      continue;
    end

    duplicateIndices = find(dt == 0);
    if any(times(duplicateIndices) ~= 0)
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        'A repeated time occurs somewhere other than t = 0.'); %#ok<AGROW>
    end
    if ~isempty(duplicateIndices)
      nDuplicateZeroTrials = nDuplicateZeroTrials + 1;
    end

    if times(end) >= totalMS
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        sprintf('Last event time %.6g ms is outside the %.6g-ms interval.', ...
        times(end), totalMS)); %#ok<AGROW>
    end

    positiveDt = dt(dt > 0);
    allPositiveIntervals = [allPositiveIntervals, positiveDt]; %#ok<AGROW>

    finalInterval = totalMS - times(end);
    finalIntervals(end+1) = finalInterval; %#ok<AGROW>

    if isfinite(row.cohNoiseFrameMS) && ...
        any(positiveDt > row.cohNoiseFrameMS + timingToleranceMS)
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'warning', ...
        sprintf('An event interval exceeds %.6g ms plus tolerance.', ...
        row.cohNoiseFrameMS)); %#ok<AGROW>
    end

    if hasPreferredNoise(k)
      if all(values == 0)
        issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
          'Preferred-noise trial contains only zero preferred-noise values.'); %#ok<AGROW>
      end
    elseif any(values ~= 0)
      issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
        ['Trial classified as having no preferred-direction noise ' ...
         'contains nonzero preferred-noise values.']); %#ok<AGROW>
    end

    observedNoiseAbs = [observedNoiseAbs, abs(values(values ~= 0))]; %#ok<AGROW>
  end

  row.nDuplicateZeroTrials = nDuplicateZeroTrials;

  if isempty(allPositiveIntervals)
    row.minPositiveEventIntervalMS = NaN;
    row.maxPositiveEventIntervalMS = NaN;
  else
    row.minPositiveEventIntervalMS = min(allPositiveIntervals);
    row.maxPositiveEventIntervalMS = max(allPositiveIntervals);
  end

  if isempty(finalIntervals)
    row.minFinalIntervalMS = NaN;
    row.maxFinalIntervalMS = NaN;
  else
    row.minFinalIntervalMS = min(finalIntervals);
    row.maxFinalIntervalMS = max(finalIntervals);
  end

  if isempty(observedNoiseAbs)
    row.maxObservedAbsNoisePC = 0;
  else
    row.maxObservedAbsNoisePC = max(observedNoiseAbs);
  end

  nonzeroPrefAmplitudes = unique(abs(prefCohNoisePC(prefCohNoisePC ~= 0)));
  row.nPreferredNoiseAmplitudes = numel(nonzeroPrefAmplitudes);
  row.preferredNoiseAmplitudesPC = ...
    string(mat2str(nonzeroPrefAmplitudes(:)'));

  if isfinite(row.headerPrefCohNoisePC) && ...
      ~isempty(nonzeroPrefAmplitudes) && ...
      any(abs(nonzeroPrefAmplitudes - abs(row.headerPrefCohNoisePC)) > 1e-6)
    issueRows(end+1) = makeIssue(fileName, NaN, 'warning', ...
      ['Trial-level preferred-noise amplitude differs from ' ...
       'sessionHeader.prefCohNoisePC.']); %#ok<AGROW>
  end

  fileIssues = strcmp({issueRows.fileName}, fileName);
  row.nErrors = sum(fileIssues & strcmp({issueRows.severity}, 'error'));
  row.nWarnings = sum(fileIssues & strcmp({issueRows.severity}, 'warning'));
  row.nIssues = row.nErrors + row.nWarnings;

  summaryRows(iFile) = row;
end

sessionTable = struct2table(summaryRows);
issueTable = struct2table(issueRows);
if ~isempty(issueTable)
  fprintf('Validation completed: %d errors and %d warnings across %d session files.\n', ...
    sum(strcmp(issueTable.severity, 'error')), sum(strcmp(issueTable.severity, 'warning')), ...
    height(sessionTable));
  disp(issueTable);
  error('dailyUpdate:BetaValidationFailed', 'validateBetaSessionData reported %d issue(s).', height(issueTable));
end

end

function row = emptySummaryRow()
row = struct( ...
  'fileName', '', ...
  'version', NaN, ...
  'nTrials', NaN, ...
  'nPreferredNoiseTrials', NaN, ...
  'nNoPreferredNoiseTrials', NaN, ...
  'nCohNoiseTrials', NaN, ...
  'nNoCohNoiseTrials', NaN, ...
  'nSignalCoherences', NaN, ...
  'signalCoherencesPC', "", ...
  'minSignalCohPC', NaN, ...
  'maxSignalCohPC', NaN, ...
  'fractionCorrect', NaN, ...
  'headerPrefCohNoisePC', NaN, ...
  'nPreferredNoiseAmplitudes', NaN, ...
  'preferredNoiseAmplitudesPC', "", ...
  'maxObservedAbsNoisePC', NaN, ...
  'preStepMS', NaN, ...
  'stepMS', NaN, ...
  'cohNoiseFrameMS', NaN, ...
  'frameRateHz', NaN, ...
  'nDuplicateZeroTrials', 0, ...
  'minPositiveEventIntervalMS', NaN, ...
  'maxPositiveEventIntervalMS', NaN, ...
  'minFinalIntervalMS', NaN, ...
  'maxFinalIntervalMS', NaN, ...
  'nIssues', 0, ...
  'nErrors', 0, ...
  'nWarnings', 0);
end

function row = emptyIssueRow()
row = struct('fileName','', 'trialIdx',NaN, 'severity','', 'message','');
end

function row = makeIssue(fileName, trialIdx, severity, message)
row = struct('fileName',fileName, 'trialIdx',trialIdx, ...
  'severity',severity, 'message',message);
end

function value = getHeaderScalar(H, fieldName)
if ~isfield(H, fieldName)
  value = NaN;
  return;
end
value = H.(fieldName);
while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
if ~(isnumeric(value) || islogical(value)) || ...
    isempty(value) || ~isscalar(value)
  value = NaN;
else
  value = double(value);
end
end
