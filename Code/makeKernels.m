  function [allProbeDirs, staleProbeDirs] = makeKernels(replace)
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
path = folderPath();
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
  clear trials;

  % Load lightweight parent metadata first.  sessionHeader is the compact
  % analysis-facing manifest used for stale-output checks.
  sessionHeader = [];
  load(dataFilePath, 'header', 'sessionHeader');

  if applyExperimentalValidityChecks(header)
    continue;
  end

  needSessionHeaderRefresh = ...
    replace || ...
    ~exist('sessionHeader', 'var') || ...
    isempty(sessionHeader) || ...
    ~isfield(sessionHeader, 'probeDirectionsDeg') || ...
    isempty(sessionHeader.probeDirectionsDeg) || ...
    ~isfield(sessionHeader, 'probeTags') || ...
    isempty(sessionHeader.probeTags);

  if needSessionHeaderRefresh
    load(dataFilePath, 'trials');
    trialMeta = trialMetaFromTrials(header, trials);
    sessionHeader = makeSessionHeader(header, trialMeta);
    save(dataFilePath, 'sessionHeader', '-append');
  end

  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;

  % ---- Tally probe directions and check whether any outputs are missing ----
  needsKernels = replace;

  for p = 1:numel(probeDirectionsDeg)
    probeDirDeg = probeDirectionsDeg(p);
    probeTag = char(probeTags{p});

    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>

    if ~replace
      analysisBaseName = sprintf('%s_%s', baseName, probeTag);

      probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
      probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
      kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
      matrixFolder = validFolder(fullfile(probeDataFolder, 'NoiseMatrices'));
      summaryFolder = validFolder(fullfile(probeDataFolder, 'KernelSummaries'));
      summaryFilePath = fullfile(summaryFolder, [analysisBaseName '_kernelSummary.mat']);

      kernelFilePath = fullfile(kernelFolder, [analysisBaseName '.mat']);
      matrixFilePath = fullfile(matrixFolder, [analysisBaseName '.mat']);
      plotFilePath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

      outputsStale = ...
        ~isfile(plotFilePath) || ...
        ~isfile(kernelFilePath) || ...
        ~isfile(matrixFilePath) || ...
        ~isfile(summaryFilePath);

      if outputsStale
        staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
        needsKernels = true;
      end
    end
  end

  if ~needsKernels
    continue
  end
  probeSessions = splitTrialsByProbeDirection(header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeDirDeg = probeSessions(p).probeDirDeg;
    probeTag = probeSessions(p).probeTag;
    sessionProbeHeader = probeSessions(p).sessionProbeHeader;
    sessionHeader = probeSessions(p).sessionHeader;
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
    
    save(matrixFilePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr','prefNoiseByPatch', 'probeNoiseByPatch', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    save(kernelFilePath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
      'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    summaryFolder = validFolder(fullfile(probeDataFolder, 'KernelSummaries'));

    compileKernelSessionSummary(kernelFilePath, ...
      'noiseFile', matrixFilePath, ...
      'summaryDir', summaryFolder, ...
      'replace', true, ...
      'doBootstrap', false);

    titleStr = sprintf('\\bf%d° Probe Kernels %s', probeDirDeg, analysisBaseName);

    plotKernels(1, titleStr, sessionHeader, kernels(1:5,:,:,:), kVars(1:5,:,:), compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
end

% function [header, sessionHeader] = ensureSessionProbeHeaders(header, sessionHeader, dataFilePath)
% % ensureSessionProbeHeaders  Cache derived probe-session metadata in parent header.
% %
% % This is a private cache in the converted file, used only to avoid loading
% % the large trials variable during stale-output checks. Downstream analysis
% % files must save individual sessionProbeHeader records, not the parent header.
% 
% if isempty(sessionHeader)
%   sessionHeader = makeSessionHeader(header);
% end
% 
% % Keep the parent-header cache aligned with the authoritative sessionHeader
% % used for derived outputs. The cache is only for fast stale-output checks;
% % probe-specific output files save their own sessionProbeHeader records.
% header.sessionHeader = sessionHeader;
% 
% if isfield(header, 'sessionProbeHeaders') && ~isempty(header.sessionProbeHeaders)
%   save(dataFilePath, 'header', 'sessionHeader', '-append');
%   return;
% end
% 
% load(dataFilePath, 'trials');
% 
% % Retain old-format compatibility at the split boundary only.
% if ~isfield(header, 'prefDirDeg')
%   header.prefDirDeg = struct('data', trials{1}.changeDots.data.directionDeg);
% end
% 
% probeSessions = splitTrialsByProbeDirection(header, trials);
% header.sessionProbeHeaders = [probeSessions.sessionProbeHeader];
% 
% save(dataFilePath, 'header', 'sessionHeader', '-append');
% end

% function tf = matFileHasVariable(filePath, varName)
% % matFileHasVariable  True if MAT file exists and contains variable varName.
% 
% tf = false;
% 
% if ~isfile(filePath)
%   return;
% end
% 
% tic
% vars = whos('-file', filePath);
% elapsed = toc;
% fprintf('%.3fs looking for %s in %s\n', elapsed, varName, filePath);
% tf = any(strcmp({vars.name}, varName));
% if tf
%   fprintf("  ***\n");
% else
%   fprintf("000\n");
% end
% end
% 
% function v = localDataValue(x)
% if isstruct(x) && isfield(x, 'data')
%   v = x.data;
% else
%   v = x;
% end
% 
% if isnumeric(v) && ~isscalar(v)
%   v = v(1);
% end
% end

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
trialMeta.probeTags = arrayfun(@(d) sprintf('probe%d', round(d)), ...
  trialMeta.probeDirectionsDeg, 'UniformOutput', false);
trialMeta.nNoiseTrials = sum(trialHasNoise & trialProbeDirs ~= -1);
end