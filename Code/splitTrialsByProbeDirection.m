function probeSessions = splitTrialsByProbeDirection(header, trials)
% splitTrialsByProbeDirection  Split one recording session into probe-specific analysis sessions.
%
% Each output element has:
%   .probeDirDeg
%   .probeTag
%   .sessionProbeHeader
%   .trials
%   .trialIdx
%
% Old single-probe files remain compatible only at this split boundary.
% If per-trial probeDirDeg is absent, header.probeDirDeg.data may be used
% to construct sessionProbeHeader. No downstream derived file should save or
% depend on the parent header.

nTrials = numel(trials);
trialProbeDirs = nan(1, nTrials);
trialHasNoise = false(1, nTrials);

% ---- Preferred new format: selected probe direction stored per trial ----
for t = 1:nTrials
  if isfield(trials{t}, 'trial') && isfield(trials{t}.trial, 'data')
    D = trials{t}.trial.data;

    if isfield(D, 'cohNoise')
      trialHasNoise(t) = logical(D.cohNoise);
    else
      % Old files did not necessarily mark no-noise trials explicitly.
      % Treat them as noise trials unless told otherwise.
      trialHasNoise(t) = true;
    end

    if isfield(D, 'probeDirDeg')
      trialProbeDirs(t) = double(D.probeDirDeg);
    end
  end
end

% ---- Backward compatibility: old single-probe sessions ----
if all(isnan(trialProbeDirs(trialHasNoise)))
  if isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
    probeDirDeg = double(header.probeDirDeg.data);
    trialProbeDirs(:) = probeDirDeg;
  else
    error('splitTrialsByProbeDirection:MissingProbeDir', ...
      ['No per-trial trials{t}.trial.data.probeDirDeg found, and no ' ...
       'header.probeDirDeg.data field is available for old-format compatibility.']);
  end
end

% ---- Assertions: fail early if partly encoded or malformed ----
if any(isnan(trialProbeDirs(trialHasNoise)))
  badNoise = find(trialHasNoise & isnan(trialProbeDirs), 1, 'first');
  error('splitTrialsByProbeDirection:IncompleteProbeDir', ...
    'Missing trial.data.probeDirDeg for noise trial %d of %d.', badNoise, nTrials);
end

if any(~isfinite(trialProbeDirs(trialHasNoise)))
  error('splitTrialsByProbeDirection:BadProbeDir', ...
    'Non-finite probe direction values found on noise trials.');
end

% Exclude no-noise trials from derived probe sessions.  cohNoise is the
% authoritative flag for whether a trial contains coherence noise.  The
% older probeDirDeg == -1 sentinel is also excluded for safety.
probeDirs = unique(trialProbeDirs(trialHasNoise));
probeDirs(probeDirs == -1) = [];
if isempty(probeDirs)
  error('splitTrialsByProbeDirection:NoProbeNoiseTrials', ...
    'No valid probe directions found after excluding probeDirDeg == -1 no-noise trials.');
end

probeSessions = repmat(struct( ...
  'probeDirDeg', [], ...
  'probeTag', '', ...
  'sessionProbeHeader', [], ...
  'trials', [], ...
  'trialIdx', []), 1, numel(probeDirs));

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialHasNoise & trialProbeDirs == probeDirDeg);
  probeTag = sprintf('probe%d', round(probeDirDeg));
  probeTrials = trials(idx);

  sessionProbeHeader = makeSessionProbeHeader(header, probeTrials, probeDirDeg, probeTag, nTrials, idx);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = probeTag;
  probeSessions(p).sessionProbeHeader = sessionProbeHeader;
  probeSessions(p).trials = probeTrials;
  probeSessions(p).trialIdx = idx;
end
end

%% makeSessionProbeHeader  Build authoritative metadata for one derived probe session.
function H = makeSessionProbeHeader(parentHeader, probeTrials, probeDirDeg, probeTag, parentNTrials, trialIdx)
%
% H is the metadata record for a probe-specific analysis session.
% Downstream kernel/noise/summary files should save H as sessionProbeHeader,
% not the original parent file header.
%
% The parent header may be used here only for file-level constants and
% provenance. Trial-varying probe metadata must be derived from probeTrials.

H = struct();

% ---- Identity / provenance ----
H.probeDirDeg = probeDirDeg;
H.probeTag = probeTag;
H.parentNumberOfTrials = parentNTrials;
H.trialIdx = trialIdx(:)';

if isfield(parentHeader, 'fileName')
  H.parentFileName = parentHeader.fileName;
end

% ---- Preferred coherence noise amplitude ----
% prefCohNoisePC is file-level by experimental contract, but we derive it
% from trials if it is absent from the parent header.  The derived value is
% stored in sessionProbeHeader, not assumed downstream from parent header.

if isfield(parentHeader, 'prefCohNoisePC')
  H.prefCohNoisePC = parentHeader.prefCohNoisePC;
elseif isfield(parentHeader, 'prefNoiseCohPC')
  H.prefCohNoisePC = parentHeader.prefNoiseCohPC;
else
  prefCohs = nan(1, numel(probeTrials));

  for t = 1:numel(probeTrials)
    D = probeTrials{t}.trial.data;

    assert(isfield(D, 'prefCohNoisePC'), ...
      'makeSessionProbeHeader:MissingPrefCohNoisePC', ...
      'Missing trial.data.prefCohNoisePC for probe trial %d.', t);

    prefCohs(t) = double(D.prefCohNoisePC);
  end

  prefCohs = prefCohs(isfinite(prefCohs));
  prefCohs = unique(round(prefCohs, 6));

  assert(numel(prefCohs) == 1, ...
    'makeSessionProbeHeader:MixedPrefCohNoisePC', ...
    'Expected exactly one prefCohNoisePC value in derived probe session.');

  H.prefCohNoisePC = struct('data', prefCohs);
end

% ---- File-level constants required downstream ----
% These are copied from the parent header because they are fixed for the
% parent data file and define the time base / kernel length.
requiredParentFields = {'frameRateHz', 'preStepMS', 'stepMS'};

for f = 1:numel(requiredParentFields)
  fieldName = requiredParentFields{f};
  assert(isfield(parentHeader, fieldName), ...
    'makeSessionProbeHeader:MissingParentField', ...
    'parentHeader.%s is required for sessionProbeHeader.', fieldName);
  H.(fieldName) = parentHeader.(fieldName);
end
% Useful if present, but do not require for all historical files.
optionalParentFields = {'prefDirDeg', 'subject', 'date', 'taskName', 'fileName'};

for f = 1:numel(optionalParentFields)
  fieldName = optionalParentFields{f};
  if isfield(parentHeader, fieldName)
    H.(fieldName) = parentHeader.(fieldName);
  end
end

% ---- Trial-derived probe direction validation ----
trialProbeDirs = nan(1, numel(probeTrials));

for t = 1:numel(probeTrials)
  D = probeTrials{t}.trial.data;
  if ~isfield(D, 'probeDirDeg')
    assert(isfield(parentHeader, 'probeDirDeg'), ...
      'makeSessionProbeHeader:MissingProbeDirDeg', ...
      'Missing parentHeader.data.probeDirDeg');
    trialProbeDirs(t) = parentHeader.probeDirDeg.data;
  else
    assert(isfield(D, 'probeDirDeg'), ...
      'makeSessionProbeHeader:MissingProbeDirDeg', ...
      'Missing trial.data.probeDirDeg for probe trial %d.', t);
    trialProbeDirs(t) = double(D.probeDirDeg);
  end
end

trialProbeDirs = unique(trialProbeDirs);

assert(numel(trialProbeDirs) == 1 && trialProbeDirs == probeDirDeg, ...
  'makeSessionProbeHeader:MixedProbeDirDeg', ...
  'Derived probe session contains inconsistent probeDirDeg values.');

% ---- Trial-derived probe coherence noise amplitude ----
probeCohs = nan(1, numel(probeTrials));

for t = 1:numel(probeTrials)
  D = probeTrials{t}.trial.data;

  assert(isfield(D, 'probeCohNoisePC'), ...
    'makeSessionProbeHeader:MissingProbeCohNoisePC', ...
    'Missing trial.data.probeCohNoisePC for probe trial %d.', t);

  probeCohs(t) = double(D.probeCohNoisePC);
end

probeCohs = probeCohs(isfinite(probeCohs));
probeCohs = unique(round(probeCohs, 6));

assert(numel(probeCohs) == 1, ...
  'makeSessionProbeHeader:MixedProbeCohNoisePC', ...
  'Expected exactly one probeCohNoisePC value in derived probe session.');

% Preserve the parameter-field style used elsewhere in the pipeline.
H.probeCohNoisePC = struct('data', probeCohs);

% ---- Noise-trial counts ----
cohNoiseFlags = true(1, numel(probeTrials));

for t = 1:numel(probeTrials)
  D = probeTrials{t}.trial.data;
  if isfield(D, 'cohNoise')
    cohNoiseFlags(t) = logical(D.cohNoise);
  end
end

H.nTrials = numel(probeTrials);
H.nNoiseTrials = sum(cohNoiseFlags);
H.nNoNoiseTrials = sum(~cohNoiseFlags);

assert(H.nNoiseTrials == H.nTrials, ...
  'makeSessionProbeHeader:UnexpectedNoNoiseTrials', ...
  'Derived probe session unexpectedly contains no-noise trials.');
end