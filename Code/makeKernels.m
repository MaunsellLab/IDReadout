  function [allProbeDirs, staleProbeDirs] = makeKernels(replace)
% makeKernels  Build and plot probe-session kernels
%
% Derived per-probe files written by makeKernels contain both headers:
%   Data/ProbeSessions/ProbeXX/Kernels
%

if nargin < 1 || isempty(replace)
  replace = false;
end
staleProbeDirs = [];
allProbeDirs = [];

% ---- Find all relevant .mat data files ----
sessionDataFolder = fullfile(folderPath(), 'Data', 'Sessions');
allMatFiles = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(allMatFiles)
  fprintf('No .mat files found in %s\n', sessionDataFolder);
  return;
end
names = {allMatFiles.name};
isFileInfo = endsWith(names, '_fileInfo.mat');
dataFiles = allMatFiles(~isFileInfo);
if isempty(dataFiles)
  fprintf('No data session files found in %s\n', sessionDataFolder);
  return;
end

% ---- Process each session file to see if probeSession kernels are missing ----
for k = 1:numel(dataFiles)
  dataFileName = dataFiles(k).name;
  dataFilePath = fullfile(sessionDataFolder, dataFileName);
  [~, baseName] = fileparts(dataFileName);

  % get a list of all probe directions using the headers in the data file
  load(dataFilePath, 'header', 'sessionHeader');
  probeDirectionsDeg = sessionHeader.probeDirectionsDeg;
  probeTags = sessionHeader.probeTags;
  needsKernels = replace;
  for p = 1:numel(probeDirectionsDeg)
    probeDirDeg = probeDirectionsDeg(p);
    allProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
    probeTag = probeTags{p};
    [kernelPath, plotPath] = makeFilePaths(probeTag, baseName);
    if ~replace
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
  % --- otherwise, load the full data set and make/plot kernels for each
  load(dataFilePath, 'trials');
  probeSessions = splitTrialsByProbeDirection(header, trials, sessionHeader);
  for p = 1:numel(probeSessions)
    probeTag = probeTags{p};
    fprintf('      processing %s [%s] ...\n', dataFileName, probeTag);
    probeDataPath = fullfile(folderPath(), 'Data', 'ProbeSessions', probeTag, ...
                      sprintf('%s_%s.mat', baseName, probeTag));

    load(probeDataPath, 'sessionProbeHeader', 'sideTypeNames', 'lr', ...
          'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', 'changeSidesAll', 'changeIndicesAll');

    % [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
    %               extractPatchNoiseMatrices(sessionHeader, sessionProbeHeader, [1 2]);

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
    
    [kernelPath, plotPath] = makeFilePaths(probeTag, baseName);
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

function [kernelPath, plotPath] = makeFilePaths(probeTag, baseName)

path = folderPath();

analysisBaseName = sprintf('%s_%s', baseName, probeTag);
plotRoot = validFolder(fullfile(path, 'Plots', 'Kernels'));
probePlotFolder = validFolder(fullfile(plotRoot, probeTag));
plotPath = fullfile(probePlotFolder, sprintf('%s.pdf', analysisBaseName));

probeDataFolder = validFolder(fullfile(path, 'Data', probeTag));
kernelFolder = validFolder(fullfile(probeDataFolder, 'Kernels'));
kernelPath = fullfile(kernelFolder, [analysisBaseName '.mat']);

end