function [sessionTable, issueTable] = validateBetaSessionData()
% validateBetaSessionData  Validate expanded BetaAnalysis SessionData files.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%
% Supports sessionNoise.version >= 2, containing all valid increment trials,
% including both preferred-noise and no-preferred-noise trials.

cleanupObj = initProjectPath(); %#ok<NASGU>

dataFolder = fullfile(folderPath(), 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');

files = dir(fullfile(dataFolder, '*.mat'));
if isempty(files)
  error('validateBetaSessionData:NoFiles', ...
    'No SessionData files found in %s.', dataFolder);
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

  requiredFields = { ...
    'version', 'sessionHeader', 'nTrials', 'trialIdx', ...
    'trialOutcome', 'changeSide', 'baseCohPC', 'stepCohPC', ...
    'signalCohPC', 'hasCohNoise', 'prefCohNoisePC', ...
    'hasPreferredNoise', 'noiseTimesMS', 'noiseCohsPC'};

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

    % if numel(times) ~= numel(values)
    %   issueRows(end+1) = makeIssue(fileName, parentTrial, 'error', ...
    %     sprintf('%d times but %d preferred-noise values.', ...
    %     numel(times), numel(values))); %#ok<AGROW>
    %   continue;
    % end

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

    % if isfinite(row.cohNoiseFrameMS) && ...
    %     (finalInterval <= 0 || ...
    %      finalInterval > row.cohNoiseFrameMS + timingToleranceMS)
    %   issueRows(end+1) = makeIssue(fileName, parentTrial, 'warning', ...
    %     sprintf(['Final represented interval is %.6g ms; expected ' ...
    %              '>0 and no longer than one noise frame.'], ...
    %     finalInterval)); %#ok<AGROW>
    % end

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

disp(sessionTable);

if isempty(issueTable)
  fprintf('Validation passed: no issues in %d session files.\n', ...
    height(sessionTable));
else
  fprintf(['Validation completed: %d errors and %d warnings ' ...
           'across %d session files.\n'], ...
    sum(strcmp(issueTable.severity, 'error')), ...
    sum(strcmp(issueTable.severity, 'warning')), ...
    height(sessionTable));
  disp(issueTable);
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
