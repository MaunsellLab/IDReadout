  function [allProbeDirs, staleProbeDirs] = makeKernels(varargin)
% makeKernels  Build and plot probe-session kernels
%
% Derived per-probe files written by makeKernels contain both headers:
%   Data/ProbeXX/Kernels
%

p = inputParser;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) isempty(x) || islogical(x));
parse(p, varargin{:});
opts = p.Results;
staleProbeDirs = [];
allProbeDirs = [];

% ---- Find all relevant session data files ----
sessionDataFolder = char(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));
dataFilePaths = selectAnalysisFiles(sessionDataFolder, 'Animal', opts.Animal);
if isempty(dataFilePaths)
  fprintf('No session files found\n');
  return;
end
% [~, sideTypeNames] = sideTypeIndex();


% ---- Find all relevant .mat data files ----
% sessionDataFolder = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions');
% allMatFiles = dir(fullfile(sessionDataFolder, '*.mat'));
% if isempty(allMatFiles)
%   fprintf('No .mat files found in %s\n', sessionDataFolder);
%   return;
% end
% names = {allMatFiles.name};
% isFileInfo = endsWith(names, '_fileInfo.mat');
% dataFiles = allMatFiles(~isFileInfo);
% if isempty(dataFiles)
%   fprintf('No data session files found in %s\n', sessionDataFolder);
%   return;
% end

% ---- Process each session file to see if probeSession kernels are missing ----
for k = 1:numel(dataFilePaths)
  dataFilePath = dataFilePaths{k};
  [~, baseName, extension] = fileparts(dataFilePath);
  dataFileName = [baseName, extension];

  % get a list of all probe directions using the headers in the data file
  load(dataFilePath, 'sessionHeader');
  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;
  needsKernels = opts.Replace;
  for p = 1:numel(probeDirectionsDeg)
    probeDirDeg = probeDirectionsDeg(p);
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
    probeTag = probeTags{p};
    [kernelPath, plotPath] = makeFilePaths(probeTag, baseName);
    if ~opts.Replace
      outputsStale = ~isfile(plotPath) || ~isfile(kernelPath);
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

  % --- otherwise, load each existing probeSession and make/plot its kernel
  for p = 1:numel(probeTags)
    probeTag = probeTags{p};
    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);
    probeDataPath = fullfile(domainFolder(mfilename('fullpath')), 'Data', probeTag, 'ProbeSessions', ...
                      sprintf('%s_%s.mat', baseName, probeTag));
    if ~isfile(probeDataPath)
      error('makeKernels:MissingProbeSession', ['Missing probeSession file:\n%s\n' ...
        'Run makeProbeSessions before makeKernels.'], probeDataPath);
    end
    load(probeDataPath, 'sessionProbeHeader', 'lr', 'prefNoiseByPatch', 'probeNoiseByPatch', ...
                'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll');
    [~, sideTypeNames] = sideTypeIndex();
    sessionData = struct;
    sessionData.sessionHeader = sessionHeader;
    sessionData.sessionProbeHeader = sessionProbeHeader;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.chosenSidesAll = chosenSidesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);
    
    [kernelPath, plotPath] = makeFilePaths(probeTag, baseName);
    save(kernelPath, 'sessionHeader', 'sessionProbeHeader', 'sideTypeNames', 'lr', 'kernels', 'kVars', 'kStats', ...
      'trialOutcomesAll', 'changeSidesAll', 'chosenSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    
    probeDirDeg = sessionProbeHeader.probeDirDeg;
    titleStr = sprintf('\\bf%d° Probe Kernels %s-%s', probeDirDeg, baseName, probeTag);
    sideTypes = [1, 2, 3, 8, 9];
    plotKernels(1, titleStr, sideTypes, sessionHeader, kernels(:,:,:,:), kVars(:,:,:), ...
        compStats, hitStats, probeDirDeg);
    exportgraphics(gcf, plotPath, 'ContentType', 'vector');
  end
end
allProbeDirs = unique(allProbeDirs);
staleProbeDirs = unique(staleProbeDirs);
  end

function [kernelPath, plotPath] = makeFilePaths(probeTag, baseName)

path = domainFolder(mfilename('fullpath'));

analysisBaseName = sprintf('%s_%s', baseName, probeTag);
plotRoot = validFolder(fullfile(path, 'Plots', 'Probes'));
probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
plotPath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
kernelPath = fullfile(kernelFolder, [analysisBaseName '.mat']);

end