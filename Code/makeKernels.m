function [allProbeDirs, staleProbeDirs] = makeKernels(replace, path)
% makeKernels  Generate kernel, noise-matrix, and plot files for converted session files.
%
%   staleProbeDirs = makeKernels(replace)
%       If replace is false (default), only computes missing outputs.
%       If replace is true, recomputes all kernels and plots even if they exist.
%
%   Directory conventions:
%       Converted data files:   path/Data/Converted/*.mat   (excluding *_fileInfo.mat)
%       Kernel .mat files:      path/Data/Kernels/<baseName>.mat
%       Noise matrix .mat files:path/Data/NoiseMatrices/<baseName>.mat
%       Kernel plot PDFs:       path/Plots/Kernels/<baseName>/<baseName>_probeXX.pdf
%
%   Output:
%       staleProbeDirs          Unique probeDirDeg values for sessions whose
%                               outputs were newly created or updated.

% ---- Handle inputs ----
if nargin < 1 || isempty(replace)
  replace = false;
end
if nargin < 2 || isempty(path)
  path = folderPath();
end
path = char(path);                % older MATLABs are happier with char for fullfile
staleProbeDirs = [];
allProbeDirs = [];

% ---- Define directories and create as needed ----
[dataFolder, existed] = validFolder(fullfile(path, 'Data', 'Converted'));
if ~existed
  error('makeKernels: MissingConvertedFolder -- Converted data folder not found: %s', dataFolder);
end
[plotRoot] = validFolder(fullfile(path, 'Plots', 'Kernels'));

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
sideTypeNames = {'diff', 'change', 'noChange', 'L', 'R', 'RF', 'Opp'};

% ---- Process each data file ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  dataFilePath = fullfile(dataFolder, dataFileName);
  [~, baseName] = fileparts(dataFileName);

  % check whether any kernels need to be made
  load(dataFilePath, 'header');
  % Ensure all the required fields are present
  header = addMissingHeaderFields(header, dataFilePath);

  % tally probe directions and check whether any needs processing
  needsKernels = replace;
  for p = 1:numel(header.probeDirsDeg)
    probeDirDeg = header.probeDirsDeg(p);
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW> % Record all probe directions represented in the data

    if ~replace
      probeTag = char(header.probeTags(p));
      analysisBaseName = sprintf('%s_%s', baseName, probeTag);

      probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
      probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
      kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
      matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));

      kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
      matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
      plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

      if ~isfile(plotFilePath) || ~isfile(kernelFilePath) || ~isfile(matrixFilePath)
        staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
        needsKernels = true;
      end
    end
  end
  if ~needsKernels
    continue
  end

  % loading trials is very slow, so we don't do it unless we need to
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbeDirection(header, trials);
  for p = 1:numel(probeSessions)
    probeDirDeg = probeSessions(p).probeDirDeg;
    probeTag = probeSessions(p).probeTag;
    probeHeader = probeSessions(p).header;
    probeTrials = probeSessions(p).trials;
    analysisBaseName = sprintf('%s_%s', baseName, probeTag);

    probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
    probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
    kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
    matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));

    kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
    matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
    plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

    fprintf('Processing %s [%s] ...\n', dataFileName, probeTag);
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
      extractPatchNoiseMatrices(probeHeader, probeTrials, [1 2]);

    lr = sessionLRMap(probeTrials);

    sessionData = struct;
    sessionData.header = probeHeader;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);

    header = probeHeader;
    save(matrixFilePath, 'header', 'sideTypeNames', 'lr','prefNoiseByPatch', 'probeNoiseByPatch', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');

    save(kernelFilePath, 'header', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', 'trialOutcomesAll', ...
      'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');

    plotKernels(1, analysisBaseName, header, kernels, kVars, compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end


%% Add key fields if they are missing from a mat file
function header = addMissingHeaderFields(header, dataFilePath)

% Check for prefDirDeg, probeDirsDeg and probeTags
if ~isfield(header, 'prefDirDeg') || ~isfield(header, 'probeDirsDeg') || ~isfield(header, 'probeTags')
  load(dataFilePath, 'trials');

  header.prefDirDeg = struct('data', trials{1}.changeDots.data.directionDeg);

  probeSessions = splitTrialsByProbeDirection(header, trials);
  probeDirsDeg = nan(1, numel(probeSessions));
  probeTags = cell(1, numel(probeSessions));
  for p = 1:numel(probeSessions)
    probeDirsDeg(p) = probeSessions(p).probeDirDeg;
    probeTags{p} = probeSessions(p).probeTag;
  end
  header.probeDirsDeg = probeDirsDeg;
  header.probeTags = probeTags;
  save(dataFilePath, 'header', '-append');
end
end
