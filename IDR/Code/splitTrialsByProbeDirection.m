function probeSessions = splitTrialsByProbeDirection(header, trials, parentSessionHeader)
% splitTrialsByProbeDirection  Split one recording session into probe-specific analysis sessions.

nTrials = numel(trials);
if nargin < 3 || isempty(parentSessionHeader)
  error('splitTrialsByProbeDirection:MissingSessionHeader', ...
    ['parentSessionHeader is required. Build sessionHeader upstream in ' ...
    'makeKernels and pass it explicitly.']);
end
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
  'sessionHeader', [], ...
  'kernelFile', [], ...
  'noiseFile', [], ...
  'trials', [], ...
  'trialIdx', []), 1, numel(probeDirs));

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialHasNoise & trialProbeDirs == probeDirDeg);
  probeTag = sprintf('Probe%d', round(probeDirDeg));
  probeTrials = trials(idx);

  sessionProbeHeader = makeSessionProbeHeader(header, parentSessionHeader, probeTrials, probeDirDeg, probeTag, ...
          nTrials, idx, probeDirs);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = probeTag;
  probeSessions(p).sessionProbeHeader = sessionProbeHeader;
  probeSessions(p).sessionHeader = parentSessionHeader;
  probeSessions(p).trials = probeTrials;
  probeSessions(p).trialIdx = idx;
end
end

%% makeSessionProbeHeader  Build authoritative metadata for one derived probe session.
function H = makeSessionProbeHeader(parentHeader, parentSessionHeader, probeTrials, probeDirDeg, probeTag, ...
parentNTrials, trialIdx, parentProbeDirectionsDeg)
%
% Create a sessionProbeHeader that describes details for one probe
% direction in a session

[~, baseName] = fileparts(parentHeader.fileName);
probeDataFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', probeTag);
kernelFolder = fullfile(probeDataFolder, 'Kernels');
noiseFolder = fullfile(probeDataFolder, 'NoiseMatrices');

H = struct();

% ---- Identity / provenance ----
H.version = 3;
H.sessionID = baseName;
H.parentFileName = parentHeader.fileName;
H.kernelFile = fullfile(kernelFolder, [baseName, '.mat']);
H.noiseFile = fullfile(noiseFolder, [baseName, '.mat']);
H.probeDirDeg = probeDirDeg;
H.probeTag = probeTag;
H.parentNumberOfTrials = parentNTrials;
H.parentNProbeDirections = numel(parentProbeDirectionsDeg);
H.parentProbeDirectionsDeg = parentProbeDirectionsDeg(:)';
H.trialIdx = trialIdx(:)';

% ---- Trial-derived probe direction validation ----
trialProbeDirs = nan(1, numel(probeTrials));

for t = 1:numel(probeTrials)
  D = probeTrials{t}.trial.data;
  if ~isfield(D, 'probeDirDeg')
    assert(isfield(parentHeader, 'probeDirDeg'), 'makeSessionProbeHeader:MissingProbeDirDeg', ...
      'Missing parentHeader.data.probeDirDeg');
    trialProbeDirs(t) = parentHeader.probeDirDeg.data;
  else
    assert(isfield(D, 'probeDirDeg'), 'makeSessionProbeHeader:MissingProbeDirDeg', ...
      'Missing trial.data.probeDirDeg for probe trial %d.', t);
    trialProbeDirs(t) = double(D.probeDirDeg);
  end
end

trialProbeDirs = unique(trialProbeDirs);

assert(isscalar(trialProbeDirs) && trialProbeDirs == probeDirDeg, ...
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

H.probeCohNoisePC = probeCohs;

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