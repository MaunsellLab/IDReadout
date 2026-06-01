  function [allProbeDirs, staleProbeDirs] = makeKernels(replace)
% makeKernels  Build all ordinary per-probe derived analysis products.
%
% makeKernels is the authoritative producer of analysis headers:
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
% Derived per-probe files written by makeKernels contain both headers:
%   Data/probeXX/NoiseMatrices
%   Data/probeXX/Kernels
%
% Aggregate files, such as average kernels and across-offset summaries, do
% not contain a single sessionHeader/sessionProbeHeader because they span
% multiple sessions.

if nargin < 1 || isempty(replace)
  replace = false;
end
staleProbeDirs = [];
allProbeDirs = [];

% ---- Define directories and create as needed ----
dataFolder = makeFilePaths();

% ---- Find all relevant .mat data files ----
allMatFiles = dir(fullfile(dataFolder, '*.mat'));
if isempty(allMatFiles)
  fprintf('No .mat files found in %s\n', dataFolder);
  return;
end
names = {allMatFiles.name};
isFileInfo = endsWith(names, '_fileInfo.mat');
dataFiles = allMatFiles(~isFileInfo);
if isempty(dataFiles)
  fprintf('No data .mat files (excluding *_fileInfo.mat) found in %s\n', dataFolder);
  return;
end
[~, sideTypeNames] = sideTypeIndex();

% ---- Process each data file ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  dataFilePath = fullfile(dataFolder, dataFileName);
  [~, baseName] = fileparts(dataFileName);
  clear header sessionHeader trials;
  % get a list of all probe directions using the headers in the data file
  load(dataFilePath, 'header', 'sessionHeader');
  needSessionHeaderRefresh = replace || ...
    ~exist('sessionHeader', 'var') || isempty(sessionHeader) || ...
    ~isfield(sessionHeader, 'probeDirectionsDeg') || isempty(sessionHeader.probeDirectionsDeg) || ...
    ~isfield(sessionHeader, 'probeTags') || isempty(sessionHeader.probeTags);
  if needSessionHeaderRefresh
    load(dataFilePath, 'trials');
    trialMeta = trialMetaFromTrials(header, trials);
    sessionHeader = makeSessionHeader(header, trialMeta);
    save(dataFilePath, 'sessionHeader', '-append');
  end
  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;

  % ---- second, tally probe directions and check whether any outputs are missing ----
  needsKernels = replace;
  for p = 1:numel(probeDirectionsDeg)
    probeDirDeg = probeDirectionsDeg(p);
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
    probeTag = char(probeTags{p});
    [~, kernelPath, noisePath, plotPath] = makeFilePaths(probeTag, baseName);
    if ~replace
      outputsStale = ~isfile(plotPath) || ~isfile(kernelPath) || ~isfile(noisePath);
      if outputsStale
        staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
        needsKernels = true;
      end
    end
  end
  % --- if nothing is missing, check the next file
  if ~needsKernels
    continue
  end
  % --- otherwise, load the full data set, split into probe dirs and make
  % kernel and noise files for each
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbeDirection(header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeTag = probeSessions(p).probeTag;
    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);

    [~, kernelPath, noisePath, plotPath] = makeFilePaths(probeTag, baseName);
    probeSessions(p).sessionProbeHeader.kernelFile = kernelPath;
    probeSessions(p).sessionProbeHeader.noiseFile = noisePath;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionHeader = probeSessions(p).sessionHeader;
    probeTrials = probeSessions(p).trials;

    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
                  extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, probeTrials, [1 2]);
    lr = sessionLRMap(probeTrials);

    sessionData = struct;
    sessionData.sessionHeader = sessionHeader;
    sessionData.sessionProbeHeader = sessionProbeHeader;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
    
    save(noisePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr','prefNoiseByPatch', 'probeNoiseByPatch', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    save(kernelPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    probeDirDeg = probeSessions(p).probeDirDeg;
    titleStr = sprintf('\\bf%d° Probe Kernels %s-%s', probeDirDeg, baseName, probeTag);
    plotKernels(1, titleStr, sessionHeader, kernels(1:5,:,:,:), kVars(1:5,:,:), compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotPath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end

function trialMeta = trialMetaFromTrials(header, trials)

nTrials = numel(trials);
trialProbeDirs = nan(1, nTrials);
trialHasNoise = false(1, nTrials);

for t = 1:nTrials
  if isfield(trials{t}, 'trial') && isfield(trials{t}.trial, 'data')
    D = trials{t}.trial.data;

    if isfield(D, 'cohNoise')
      trialHasNoise(t) = logical(D.cohNoise);
    else
      trialHasNoise(t) = true;
    end

    if isfield(D, 'probeDirDeg')
      trialProbeDirs(t) = double(D.probeDirDeg);
    end
  end
end

% Old single-probe compatibility at the conversion boundary only.
if all(isnan(trialProbeDirs(trialHasNoise)))
  if isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
    trialProbeDirs(:) = double(header.probeDirDeg.data);
  else
    error('trialMetaFromTrials:MissingProbeDir', ...
      'No per-trial probeDirDeg and no header.probeDirDeg.data.');
  end
end

if any(isnan(trialProbeDirs(trialHasNoise)))
  error('trialMetaFromTrials:IncompleteProbeDir', ...
    'Missing probeDirDeg for at least one noise trial.');
end

probeDirectionsDeg = unique(trialProbeDirs(trialHasNoise));
probeDirectionsDeg(probeDirectionsDeg == -1) = [];

trialMeta = struct();
trialMeta.probeDirectionsDeg = probeDirectionsDeg(:)';
trialMeta.nProbeDirections = numel(probeDirectionsDeg);
trialMeta.probeTags = arrayfun(@(d) sprintf('Probe%d', round(d)), trialMeta.probeDirectionsDeg, 'UniformOutput', false);
trialMeta.nNoiseTrials = sum(trialHasNoise & trialProbeDirs ~= -1);
end

function [dataFolder, kernelPath, noisePath, plotPath] = makeFilePaths(probeTag, baseName)

path = folderPath();
dataFolder = validFolder(fullfile(path, 'Data', 'Sessions'));

if nargin < 2
  kernelPath = [];
  noisePath = [];
  plotPath = [];
  return;           % only process dataFolder and plotRoot
end

analysisBaseName = sprintf('%s_%s', baseName, probeTag);

plotRoot = validFolder(fullfile(path, 'Plots', 'Kernels'));
probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
plotPath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));
kernelPath = fullfile(kernelFolder, [analysisBaseName '.mat']);
noisePath = fullfile(matrixFolder, [analysisBaseName '.mat']);

end