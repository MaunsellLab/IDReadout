function [allProbeDirs, staleProbeDirs] = makeProbeSessions(replace)
% makeProbeSessions  Split sessions into separate instances for each probe offset.
%
% makeProbeSessions is the authoritative producer of analysis headers:
%
%   sessionHeader
%       Parent-session metadata shared by all probe splits.
%       Includes timing, preferred direction/noise amplitude, and the
%       session-level probe manifest:
%           probeDirectionsDeg
%           probeTags
%
%   sessionProbeHeader
%       Probe-specific metadata for one derived probe session.
%       Includes probeDirDeg, probeTag, probeCohNoisePC, and probe-session trial counts.
%
% Derived per-probe files written by makeProbeSessions contain both headers:
%   Data/ProbeXX/ProbeSessions
cleanupObj = initProjectPath(); %#ok<NASGU>
if nargin < 1 || isempty(replace)
  replace = false;
end
staleProbeDirs = [];
allProbeDirs = [];

% ---- Find all relevant session data files ----
sessionDataFolder = fullfile(folderPath(), 'Data', 'FullSessions');
allMatFiles = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(allMatFiles)
  fprintf('No session files found\n');
  return;
end
names = {allMatFiles.name};
isFileInfo = endsWith(names, '_fileInfo.mat');
dataFiles = allMatFiles(~isFileInfo);
if isempty(dataFiles)
  fprintf('No session files found\n');
  return;
end
[~, sideTypeNames] = sideTypeIndex();

% ---- Process each data file ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  [~, baseName] = fileparts(dataFileName);
  dataFilePath = fullfile(sessionDataFolder, dataFileName);
  
  % get a list of all probe directions using the headers in the data file.
  % Don't load trials because it is too slow when we are examining whether
  % processing is needed
  clear header sessionHeader
  load(dataFilePath, 'header', 'sessionHeader');
  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;
  needsProbeSessions = replace;
  % ---- check whether this file is missing probeSessions ----
  if ~needsProbeSessions
    for p = 1:numel(probeDirectionsDeg)
      probeDirDeg = probeDirectionsDeg(p);
      allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
      probeTag = probeTags{p};
      probeDataFolder = validFolder(fullfile(folderPath(), 'Data', probeTag, 'ProbeSessions'));
      probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
      if isfile(probeSessionPath) && ~replace
        continue;
      end
      staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
      needsProbeSessions = true;
    end
  end
  % --- if nothing is missing, check the next file
  if ~needsProbeSessions
    continue
  end

  %--- otherwise, load the full data set, make probeSession files
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbe(header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeTag = probeSessions(p).probeTag;
    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);
    probeDataFolder = validFolder(fullfile(folderPath(), 'Data', probeTag, 'ProbeSessions'));
    probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
    if isfile(probeSessionPath) && ~replace
      continue;
    end
    sessionHeader = probeSessions(p).sessionHeader;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionProbeHeader.probeSessionPath = probeSessionPath;
    probeTrials = probeSessions(p).probeTrials;
    lr = sessionLRMap(probeTrials);
    trialIdx = probeSessions(p).trialIdx';
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
                  extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, probeTrials, [1 2]);

    save(probeSessionPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'probeTrials', 'trialIdx', 'lr', ...
      'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', '-v7.3');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end

%% splitTrialsByProbe
function probeSessions = splitTrialsByProbe(header, trials, parentSessionHeader)
% splitTrialsByProbeDirection  Split one recording session into probe-specific analysis sessions.

nTrials = numel(trials);
trialProbeDirs = nan(1, nTrials);
trialHasNoise = false(1, nTrials);

% ---- selected probe direction stored per trial ----
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
  end
end

% ---- Assertions: fail early if partly encoded or malformed ----
if any(isnan(trialProbeDirs(trialHasNoise)))
  badNoise = find(trialHasNoise & isnan(trialProbeDirs), 1, 'first');
  error('splitTrialsByProbe:IncompleteProbeDir', ...
    'Missing trial.data.probeDirDeg for noise trial %d of %d.', badNoise, nTrials);
end
if any(~isfinite(trialProbeDirs(trialHasNoise)))
  error('splitTrialsByProbe:BadProbeDir', ...
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
  'probeSessionPath', [], ...
  'probeTrials', [], ...
  'trialIdx', []), 1, numel(probeDirs));
% 'trials', [], ...

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialHasNoise & trialProbeDirs == probeDirDeg);
  probeTag = sprintf('Probe%d', round(probeDirDeg));
  probeTrials = trials(idx);

  sessionProbeHeader = makeSessionProbeHeader(header, probeTrials, probeDirDeg, probeTag, nTrials, probeDirs);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = probeTag;
  probeSessions(p).sessionProbeHeader = sessionProbeHeader;
  probeSessions(p).sessionHeader = parentSessionHeader;
  probeSessions(p).probeTrials = probeTrials;
  probeSessions(p).trialIdx = idx;
end
end

%% makeSessionProbeHeader  Build authoritative metadata for one derived probe session.
function H = makeSessionProbeHeader(parentHeader, probeTrials, probeDirDeg, probeTag, ...
parentNTrials, parentProbeDirectionsDeg)
%
% Create a sessionProbeHeader that describes details for one probe
% direction in a session

[~, baseName] = fileparts(parentHeader.fileName);
probeDataFolder = fullfile(folderPath, 'Data', probeTag);
probeSessionFolder = fullfile(probeDataFolder, 'ProbeSessions');

H = struct();

% ---- Identity / provenance ----
H.version = 3;
H.sessionID = baseName;
H.parentFileName = parentHeader.fileName;
H.probeSessionPath = fullfile(probeSessionFolder, [baseName, '.mat']);
H.probeDirDeg = probeDirDeg;
H.probeTag = probeTag;
H.parentNumberOfTrials = parentNTrials;
H.parentNProbeDirections = numel(parentProbeDirectionsDeg);
H.parentProbeDirectionsDeg = parentProbeDirectionsDeg(:)';

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

assert(isscalar(probeCohs), ...
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