function staleProbeDirs = makeProbeSessions(varargin)
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

p = inputParser;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) isempty(x) || islogical(x));
parse(p, varargin{:});
opts = p.Results;
staleProbeDirs = [];

% ---- Find all relevant session data files ----
sessionDataFolder = char(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));
dataFilePaths = selectAnalysisFiles(sessionDataFolder, 'Animal', opts.Animal);
if isempty(dataFilePaths)
  fprintf('No session files found\n');
  return;
end
[~, sideTypeNames] = sideTypeIndex();

% ---- Process each data file ----
for k = 1:numel(dataFilePaths)
  dataFilePath = dataFilePaths{k};
  [~, baseName] = fileparts(dataFilePath);
  
  % get a list of all probe directions using the headers in the data file.
  % Don't load trials because it is too slow when we are examining whether
  % processing is needed
  clear header sessionHeader
  load(dataFilePath, 'header', 'sessionHeader');
  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;
  needsProbeSessions = opts.Replace;
  % if we're not doing all, check whether this file is missing probeSessions
  if ~needsProbeSessions
    for p = 1:numel(probeDirectionsDeg)
      probeDirDeg = probeDirectionsDeg(p);
      probeTag = probeTags{p};
      probeDataFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', probeTag, 'ProbeSessions'));
      probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
      if isfile(probeSessionPath)
        continue;
      end
      needsProbeSessions = true;
    end
  end
  % --- if nothing is missing, check the next file
  if ~needsProbeSessions
    continue
  end

  %--- otherwise, load the full data set, make probeSession files
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbe(baseName, header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeTag = probeSessions(p).probeTag;
    fprintf('      processing %s [%s] ...\n', baseName, probeTag);
    probeDataFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', probeTag, 'ProbeSessions'));
    probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
    if isfile(probeSessionPath) && ~opts.Replace
      continue;
    end
    sessionHeader = probeSessions(p).sessionHeader;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionProbeHeader.probeSessionPath = probeSessionPath;
    probeTrials = probeSessions(p).probeTrials;
    lr = sessionLRMap(probeTrials);
    trialIdx = probeSessions(p).trialIdx';
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, chosenSidesAll, changeIndicesAll] = ...
                  extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, probeTrials, [1 2]);

    save(probeSessionPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'probeTrials', 'trialIdx', 'lr', ...
      'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', ...
      'changeIndicesAll', '-v7.3');
    staleProbeDirs = [staleProbeDirs, probeSessions(p).sessionProbeHeader.probeDirDeg]; %#ok<AGROW>
  end
end
staleProbeDirs = unique(staleProbeDirs);
end

%% splitTrialsByProbe
function probeSessions = splitTrialsByProbe(baseName, header, trials, parentSessionHeader)
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

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialHasNoise & trialProbeDirs == probeDirDeg);
  probeTag = sprintf('Probe%d', round(probeDirDeg));
  probeTrials = trials(idx);

  sessionProbeHeader = makeSessionProbeHeader(baseName, header, probeTrials, probeDirDeg, probeTag, nTrials, probeDirs);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = probeTag;
  probeSessions(p).sessionProbeHeader = sessionProbeHeader;
  probeSessions(p).sessionHeader = parentSessionHeader;
  probeSessions(p).probeTrials = probeTrials;
  probeSessions(p).trialIdx = idx;
end
end

%% makeSessionProbeHeader  Build authoritative metadata for one derived probe session.
function H = makeSessionProbeHeader(matFileName, parentHeader, probeTrials, probeDirDeg, probeTag, ...
parentNTrials, parentProbeDirectionsDeg)
%
% Create a sessionProbeHeader that describes details for one probe
% direction in a session

[~, baseName] = fileparts(parentHeader.fileName);
probeDataFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', probeTag);
probeSessionFolder = fullfile(probeDataFolder, 'ProbeSessions');

H = struct();

% ---- Identity / provenance ----
H.version = 3;
H.sessionID = baseName;
H.parentFileName = matFileName;
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

  assert(isfield(D, 'probeCohNoisePC'), 'makeSessionProbeHeader:MissingProbeCohNoisePC', ...
    'Missing trial.data.probeCohNoisePC for probe trial %d.', t);

  probeCohs(t) = double(D.probeCohNoisePC);
end

probeCohs = probeCohs(isfinite(probeCohs));
probeCohs = unique(round(probeCohs, 6));

assert(isscalar(probeCohs), 'makeSessionProbeHeader:MixedProbeCohNoisePC', ...
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

assert(H.nNoiseTrials == H.nTrials, 'makeSessionProbeHeader:UnexpectedNoNoiseTrials', ...
  'Derived probe session unexpectedly contains no-noise trials.');
end

%%================================================
function [prefNoise, probeNoise, trialOutcomes, changeSides, chosenSides, changeIndices] = ...
  extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, trials, stepTypes)
% extractPatchNoiseMatrices
% Return patch-wise noise matrices for all valid trials.
%
% OUTPUTS
%   prefNoise      : 2 x m x nTrials   (1=RF, 2=Opp)
%   probeNoise     : 2 x m x nTrials   (1=RF, 2=Opp)
%                    Effective probe noise. For paired yoked probes, the
%                    single-stream probe coherence noise is multiplied by
%                    nYokedProbeStreams before storage.
%   trialOutcomes  : 1 x nTrials       (0=correct, 1=wrong)
%   changeSides    : 1 x nTrials       (0=RF changed, 1=Opp changed)
%   chosenSides    : 1 x nTrials       (0=RF chosen, 1=Opp chosen)
%   changeIndices  : 1 x nTrials       (1=DEC, 2=INC)

% trials = sessionProbeHeader.probeTrials;
nTrials = numel(trials);
if nTrials == 0
  error('extractPatchNoiseMatrices:EmptyInput', 'Input "trials" is empty.');
end

validIdx = [];
trialOutcomes = [];
changeSides = [];
chosenSides = [];
changeIndices = [];

for k = 1:nTrials
  tr = trials{k};
  if ~isfield(tr, 'trialEnd') || ~isfield(tr, 'trialCertify') || ~isfield(tr, 'trial')
    continue;
  end

  tCert = tr.trialCertify.data;
  tEnd  = tr.trialEnd.data;
  tStep = tr.trial.data.changeIndex + 1;   % 1=DEC, 2=INC
  if ~(tCert == 0 && ismember(tEnd, [0 1]) && ismember(tStep, stepTypes))
    continue;
  end

  if ~tr.trial.data.cohNoise
    % fprintf('skipping trial %d, cohNoise %d, probeDirDeg %d\n', k, tr.trial.data.cohNoise, tr.trial.data.probeDirDeg);
    continue;
  end

  validIdx(end+1)       = k; %#ok<AGROW>
  trialOutcomes(end+1)  = tEnd; %#ok<AGROW>
  changeSides(end+1)    = tr.trial.data.changeSide; %#ok<AGROW>
  if ~isfield(tr, 'targetChosen') 
    fprintf('Bad\n');
  end
  chosenSides(end+1)    = tr.targetChosen.data; %#ok<AGROW>
  changeIndices(end+1)  = tStep; %#ok<AGROW>
end

nValid = numel(validIdx);
if nValid == 0
  error('extractPatchNoiseMatrices:NoValidTrials', 'No valid matching trials were found.');
end

% frameRateHz = localDataValue(sessionHeader.frameRateHz);
frameRateHz = sessionHeader.frameRateHz;
msPerVFrame = 1000.0 / frameRateHz;
m = round((sessionHeader.preStepMS + sessionHeader.stepMS) / msPerVFrame);

nYokedProbeStreams = probeStreamCountFromSessionProbeHeader(sessionProbeHeader);

prefNoise  = nan(2, m, nValid);
probeNoise = nan(2, m, nValid);

for kk = 1:nValid
  tr = trials{validIdx(kk)};

  prefChange    = fillFromTimes(tr.changePrefCohsPC.data(:),    tr.changeTimesMS.data(:),   m, msPerVFrame);
  probeChange   = fillFromTimes(tr.changeProbeCohsPC.data(:),   tr.changeTimesMS.data(:),   m, msPerVFrame);
  prefNoChange  = fillFromTimes(tr.noChangePrefCohsPC.data(:),  tr.noChangeTimesMS.data(:), m, msPerVFrame);
  probeNoChange = fillFromTimes(tr.noChangeProbeCohsPC.data(:), tr.noChangeTimesMS.data(:), m, msPerVFrame);
  
  % The stored probe coherence stream is the scalar used for each member of
  % the yoked pair. For kernel estimation, probeNoise should represent the
  % effective perturbation delivered by the paired probe streams.
  probeChange   = nYokedProbeStreams * probeChange;
  probeNoChange = nYokedProbeStreams * probeNoChange;

  if tr.trial.data.changeSide == 0
    % RF changed, Opp noChange
    prefNoise(1,:,kk)  = prefChange;
    probeNoise(1,:,kk) = probeChange;
    prefNoise(2,:,kk)  = prefNoChange;
    probeNoise(2,:,kk) = probeNoChange;
  else
    % Opp changed, RF noChange
    prefNoise(1,:,kk)  = prefNoChange;
    probeNoise(1,:,kk) = probeNoChange;
    prefNoise(2,:,kk)  = prefChange;
    probeNoise(2,:,kk) = probeChange;
  end
end
end

function v = fillFromTimes(cohsPC, timesMS, m, msPerVFrame)
v = nan(m,1);
if isempty(timesMS) || isempty(cohsPC)
  return;
end

nTimes = min(numel(timesMS), numel(cohsPC));
timesMS = timesMS(1:nTimes);
cohsPC  = cohsPC(1:nTimes);

for tIndex = 1:nTimes
  t0 = timesMS(tIndex);
  theVFrame = floor(t0 / msPerVFrame) + 1;
  if theVFrame < 1
    theVFrame = 1;
  elseif theVFrame > m
    continue;
  end

  if tIndex < nTimes
    t1 = timesMS(tIndex + 1);
    nextVFrame = floor(t1 / msPerVFrame) + 1;
  else
    nextVFrame = m + 1;
  end

  if nextVFrame <= theVFrame
    continue;
  end
  if nextVFrame > m + 1
    nextVFrame = m + 1;
  end

  v(theVFrame:nextVFrame-1) = cohsPC(tIndex);
end
end

%%=========================================================================
function n = probeStreamCountFromSessionProbeHeader(sessionProbeHeader)
% Number of yoked probe streams represented by the effective probe-noise
% variable.
%
%   0 < probeDirDeg < 180 : paired yoked streams at +/- probeDirDeg
%   probeDirDeg == 180   : legacy single opposite-direction stream

% assert(isfield(sessionProbeHeader, 'probeDirDeg'), ...
%   'extractPatchNoiseMatrices:MissingProbeDir', ...
%   'Cannot determine probe stream count because sessionProbeHeader.probeDirDeg is missing.');

% probeDirDeg = localHeaderScalar(sessionProbeHeader.probeDirDeg);
probeDirDeg = abs(double(sessionProbeHeader.probeDirDeg));

if probeDirDeg > 0 && probeDirDeg < 180
  n = 2;
elseif abs(probeDirDeg - 180) < 1e-9
  n = 1;
else
  error('extractPatchNoiseMatrices:UnsupportedProbeDir', 'Unsupported probeDirDeg for probe noise extraction: %g.', ...
      probeDirDeg);
end
end

%%=========================================================================
% function x = localHeaderScalar(v)
% % Accept either modern scalar fields or older struct-with-data fields.
% 
% if isstruct(v) && isfield(v, 'data')
%   v = v.data;
% end
% 
% x = v(1);
% end
% 
% %%=========================================================================
% function v = localDataValue(x)
% if isstruct(x) && isfield(x, 'data')
%   v = x.data;
% else
%   v = x;
% end
% if isnumeric(v) && ~isscalar(v)
%   v = v(1);
% end
% end