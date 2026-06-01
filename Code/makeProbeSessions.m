function makeProbeSessions(replace)
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
%       Includes probeDirDeg, probeTag, trialIdx, probeCohNoisePC, and
%       probe-session trial counts.
%
% Derived per-probe files written by makeProbeSessions contain both headers:
%   Data/ProbeSessions/ProbeXX/ProbeSessions

if nargin < 1 || isempty(replace)
  replace = false;
end

% ---- Find all relevant session data files ----
sessionDataFolder = fullfile(folderPath(), 'Data', 'Sessions');
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
      probeTag = probeTags{p};
      probeDataFolder = validFolder(fullfile(folderPath(), 'Data', 'ProbeSessions', probeTag));
      probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
      if isfile(probeSessionPath) && ~replace
        continue;
      end
      needsProbeSessions = true;
      break;
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
    probeDataFolder = validFolder(fullfile(folderPath(), 'Data', 'ProbeSessions', probeTag));
    probeSessionPath = fullfile(probeDataFolder, [sprintf('%s_%s.mat', baseName, probeTag)]);
    if isfile(probeSessionPath) && ~replace
      continue;
    end
    sessionHeader = probeSessions(p).sessionHeader;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionProbeHeader.probeSessionPath = probeSessionPath;
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
                  extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, [1 2]);
    lr = sessionLRMap(sessionProbeHeader.probeTrials);

    % sessionData = struct;
    % sessionData.sessionHeader = sessionHeader;
    % sessionData.sessionProbeHeader = sessionProbeHeader;
    % sessionData.sideTypeNames = sideTypeNames;
    % sessionData.lr = lr;
    % sessionData.prefNoiseByPatch = prefNoiseByPatch;
    % sessionData.probeNoiseByPatch = probeNoiseByPatch;
    % sessionData.trialOutcomesAll = trialOutcomesAll;
    % sessionData.changeSidesAll = changeSidesAll;
    % sessionData.changeIndicesAll = changeIndicesAll;

    % [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
    
    save(probeSessionPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'prefNoiseByPatch', ...
      'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', '-v7.3');
    
    % save(kernelPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
    %   'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    % probeDirDeg = probeSessions(p).probeDirDeg;
    % titleStr = sprintf('\\bf%d° Probe Kernels %s-%s', probeDirDeg, baseName, probeTag);
    % plotKernels(1, titleStr, sessionHeader, kernels(1:5,:,:,:), kVars(1:5,:,:), compStats, hitStats, probeDirDeg);
    % exportgraphics(gcf, plotPath, 'ContentType', 'vector');
  end
end
end

% function trialMeta = trialMetaFromTrials(header, trials)
% 
% nTrials = numel(trials);
% trialProbeDirs = nan(1, nTrials);
% trialHasNoise = false(1, nTrials);
% 
% for t = 1:nTrials
%   if isfield(trials{t}, 'trial') && isfield(trials{t}.trial, 'data')
%     D = trials{t}.trial.data;
% 
%     if isfield(D, 'cohNoise')
%       trialHasNoise(t) = logical(D.cohNoise);
%     else
%       trialHasNoise(t) = true;
%     end
% 
%     if isfield(D, 'probeDirDeg')
%       trialProbeDirs(t) = double(D.probeDirDeg);
%     end
%   end
% end
% 
% % Old single-probe compatibility at the conversion boundary only.
% if all(isnan(trialProbeDirs(trialHasNoise)))
%   if isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
%     trialProbeDirs(:) = double(header.probeDirDeg.data);
%   else
%     error('trialMetaFromTrials:MissingProbeDir', ...
%       'No per-trial probeDirDeg and no header.probeDirDeg.data.');
%   end
% end
% 
% if any(isnan(trialProbeDirs(trialHasNoise)))
%   error('trialMetaFromTrials:IncompleteProbeDir', ...
%     'Missing probeDirDeg for at least one noise trial.');
% end
% 
% probeDirectionsDeg = unique(trialProbeDirs(trialHasNoise));
% probeDirectionsDeg(probeDirectionsDeg == -1) = [];
% 
% trialMeta = struct();
% trialMeta.probeDirectionsDeg = probeDirectionsDeg(:)';
% trialMeta.nProbeDirections = numel(probeDirectionsDeg);
% trialMeta.probeTags = arrayfun(@(d) sprintf('Probe%d', round(d)), trialMeta.probeDirectionsDeg, 'UniformOutput', false);
% trialMeta.nNoiseTrials = sum(trialHasNoise & trialProbeDirs ~= -1);
% end

% %% makeFilesPaths
% function [dataFolder, kernelPath, probeSessionPath, plotPath] = makeFilePaths(probeTag, baseName)
% 
% path = folderPath();
% dataFolder = validFolder(fullfile(path, 'Data', 'Sessions'));
% 
% if nargin < 2
%   kernelPath = [];
%   probeSessionPath = [];
%   plotPath = [];
%   return;           % only process dataFolder and plotRoot
% end
% 
% analysisBaseName = sprintf('%s_%s', baseName, probeTag);
% 
% plotRoot = validFolder(fullfile(path, 'Plots', 'Kernels'));
% probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
% plotPath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));
% 
% probeDataFolder = validFolder(fullfile(path, 'Data', 'ProbeSessions', probeTag));
% kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
% matrixFolder = validFolder(fullfile(probeDataFolder, 'ProbeSessions'));
% kernelPath = fullfile(kernelFolder, [analysisBaseName '.mat']);
% probeSessionPath = fullfile(matrixFolder, [analysisBaseName '.mat']);
% 
% end

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
  'probeSessionPath', [], ...
  'trials', [], ...
  'trialIdx', []), 1, numel(probeDirs));

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialHasNoise & trialProbeDirs == probeDirDeg);
  probeTag = sprintf('Probe%d', round(probeDirDeg));
  probeTrials = trials(idx);

  sessionProbeHeader = makeSessionProbeHeader(header, probeTrials, probeDirDeg, probeTag, nTrials, idx, probeDirs);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = probeTag;
  probeSessions(p).sessionProbeHeader = sessionProbeHeader;
  probeSessions(p).sessionHeader = parentSessionHeader;
  probeSessions(p).trialIdx = idx;
end
end

%% makeSessionProbeHeader  Build authoritative metadata for one derived probe session.
function H = makeSessionProbeHeader(parentHeader, probeTrials, probeDirDeg, probeTag, ...
parentNTrials, trialIdx, parentProbeDirectionsDeg)
%
% Create a sessionProbeHeader that describes details for one probe
% direction in a session

[~, baseName] = fileparts(parentHeader.fileName);
probeDataFolder = fullfile(folderPath, 'Data', probeTag);
probeSessionFolder = fullfile(probeDataFolder, 'NoiseMatrices');

H = struct();

% ---- Identity / provenance ----
H.version = 3;
H.sessionID = baseName;
H.parentFileName = parentHeader.fileName;
H.probeSessionPath = fullfile(probeSessionFolder, [baseName, '.mat']);
H.probeDirDeg = probeDirDeg;
H.probeTag = probeTag;
H.probeTrials = probeTrials;
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