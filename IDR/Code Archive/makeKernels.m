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
  path = domainFolder(mfilename('fullpath'));
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
[~, sideTypeNames] = sideTypeIndex();

% ---- Process each data file ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  dataFilePath = fullfile(dataFolder, dataFileName);
  [~, baseName] = fileparts(dataFileName);

  % check whether any kernels need to be made
  load(dataFilePath, 'header');
  % Ensure all the required fields are present
  header = ensureSessionProbeHeaders(header, dataFilePath);
  sessionProbeHeaders = header.sessionProbeHeaders;
  if excludeFile(header)
    continue;
  end

  % tally probe directions and check whether any needs processing
  needsKernels = replace;

  for p = 1:numel(sessionProbeHeaders)
    probeDirDeg = sessionProbeHeaders(p).probeDirDeg;
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
  
    if ~replace
      probeTag = char(sessionProbeHeaders(p).probeTag);
      analysisBaseName = sprintf('%s_%s', baseName, probeTag);

      probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
      probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
      kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
      matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));

      kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
      matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
      plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

      outputsStale = ...
        ~isfile(plotFilePath) || ...
        ~isfile(kernelFilePath) || ...
        ~isfile(matrixFilePath) || ...
        ~matFileHasVariable(kernelFilePath, 'sessionProbeHeader') || ...
        ~matFileHasVariable(matrixFilePath, 'sessionProbeHeader');
      
      if outputsStale
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
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    probeTrials = probeSessions(p).trials;
    analysisBaseName = sprintf('%s_%s', baseName, probeTag);

    probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
    probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
    kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
    matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));

    kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
    matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
    plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
      extractPatchNoiseMatrices(sessionProbeHeader, probeTrials, [1 2]);
    lr = sessionLRMap(probeTrials);

    sessionData = struct;
    sessionData.sessionProbeHeader = sessionProbeHeader;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
    
    save(matrixFilePath, 'sessionProbeHeader', 'sideTypeNames', 'lr','prefNoiseByPatch', 'probeNoiseByPatch', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    save(kernelFilePath, 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    titleStr = sprintf('\\bf%d° Probe Kernels %s', probeDirDeg, analysisBaseName);

    plotKernels(1, titleStr, sessionProbeHeader, kernels(1:5,:,:,:), kVars(1:5,:,:), compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end

%% Add fields
function header = ensureSessionProbeHeaders(header, dataFilePath)
% ensureSessionProbeHeaders  Cache derived probe-session metadata in parent header.
%
% This is a private cache in the converted file, used only to avoid loading
% the large trials variable during stale-output checks. Downstream analysis
% files must save individual sessionProbeHeader records, not the parent header.

if isfield(header, 'sessionProbeHeaders') && ~isempty(header.sessionProbeHeaders)
  return;
end

load(dataFilePath, 'trials');

% Retain old-format compatibility by ensuring prefDirDeg is available if
% historical downstream/sessionProbeHeader creation needs it.
if ~isfield(header, 'prefDirDeg')
  header.prefDirDeg = struct('data', trials{1}.changeDots.data.directionDeg);
end

probeSessions = splitTrialsByProbeDirection(header, trials);
header.sessionProbeHeaders = [probeSessions.sessionProbeHeader];

save(dataFilePath, 'header', '-append');
end

function tf = matFileHasVariable(filePath, varName)
% matFileHasVariable  True if MAT file exists and contains variable varName.

tf = false;

if ~isfile(filePath)
  return;
end

vars = whos('-file', filePath);
tf = any(strcmp({vars.name}, varName));
end