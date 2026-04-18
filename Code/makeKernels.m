function staleProbeDirs = makeKernels(replace, path)
% makeKernels  Generate kernel, noise-matrix, and plot files for converted session files.
%
%   staleProbeDirs = makeKernels()
%       Uses folderPath() to get the base path.
%
%   staleProbeDirs = makeKernels(replace)
%       If replace is true, recomputes kernels and plots even if they exist.
%       If replace is false (default), only computes missing outputs.
%
%   staleProbeDirs = makeKernels(replace, path)
%       Uses the specified base path.
%
%   Directory conventions:
%       Converted data files:   path/Data/Converted/*.mat   (excluding *_fileInfo.mat)
%       Kernel .mat files:      path/Data/Kernels/<baseName>.mat
%       Noise matrix .mat files:path/Data/NoiseMatrices/<baseName>.mat
%       Kernel plot PDFs:       path/Plots/Kernels/<baseName>_probeXX.pdf
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
  path = char(path);   % older MATLABs are happier with char for fullfile

  staleProbeDirs = [];

  % ---- Define directories and create as needed ----
  dataFolder   = fullfile(path, 'Data', 'Converted');
  kernelFolder = fullfile(path, 'Data', 'Kernels');
  matrixFolder = fullfile(path, 'Data', 'NoiseMatrices');
  plotRoot = fullfile(path, 'Plots', 'Kernels');

  if ~exist(dataFolder, 'dir')
    error('makeKernels: MissingConvertedFolder -- Converted data folder not found: %s', dataFolder);
  end
  if ~exist(plotRoot, 'dir')
    mkdir(plotRoot);
  end
  if ~exist(kernelFolder, 'dir')
    mkdir(kernelFolder);
  end
  if ~exist(matrixFolder, 'dir')
    mkdir(matrixFolder);
  end

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

  % ---- Process each data file ----
  sideTypeNames = {'diff', 'change', 'noChange', 'L', 'R', 'RF', 'Opp'};
  numSkipped = 0;
  for k = 1:numel(dataFiles)
    dataFileName = dataFiles(k).name;
    dataFilePath = fullfile(dataFolder, dataFileName);
    [~, baseName] = fileparts(dataFileName);

    % ---- Load header first so probeDirDeg is available for validation and naming ----

    S = load(dataFilePath, 'header');
    header = S.header;
    probeDirDeg = header.probeDirDeg.data;
    probeTag = sprintf('probe%d', round(probeDirDeg));
    probePlotFolder = fullfile(plotRoot, probeTag);
    if ~exist(probePlotFolder, 'dir')
      mkdir(probePlotFolder);
    end
    plotFilePath = fullfile(probePlotFolder, sprintf('%s_%s.pdf', baseName, probeTag));
    if ~replace && isfile(plotFilePath)
      numSkipped = numSkipped + 1;
      continue;
    end
    kernelFilePath = fullfile(kernelFolder, [baseName '.mat']);
    matrixFilePath = fullfile(matrixFolder, [baseName '.mat']);
    if ~replace && isfile(kernelFilePath) && isfile(matrixFilePath)
      numSkipped = numSkipped + 1;
      continue;
    end
    fprintf('Processing %s [probe %g] ...\n', dataFileName, probeDirDeg);
    S = load(dataFilePath, 'trials');
    trials = S.trials;

    % For a while we weren't putting prefDirDeg into the file header. It can be
    % extracted, but it is very useful to have it available in the header.
    % This was fixed in IDR eventually, but here we update the header for older
    % data files
    if ~isfield(header, 'prefDirDeg')
      header.prefDirDeg = struct('data', trials{1}.changeDots.data.directionDeg);
      fprintf('%s updating header\n', dataFilePath);
      save(dataFilePath, 'header', 'trials');
    end

    % ---- Extract patchwise noise matrices and trial labels ----
    [prefNoiseByPatch, probeNoiseByPatch, trialOutcomesAll, changeSidesAll, changeIndicesAll] = ...
      extractPatchNoiseMatrices(header, trials, [1 2]);

    % ---- Determine sessionwise L/R mapping ----
    lr = sessionLRMap(trials);

    % ---- Package session data for canonical kernel computation ----
    sessionData = struct;
    sessionData.header = header;
    sessionData.sideTypeNames = sideTypeNames;
    sessionData.lr = lr;
    sessionData.prefNoiseByPatch = prefNoiseByPatch;
    sessionData.probeNoiseByPatch = probeNoiseByPatch;
    sessionData.trialOutcomesAll = trialOutcomesAll;
    sessionData.changeSidesAll = changeSidesAll;
    sessionData.changeIndicesAll = changeIndicesAll;

    % ---- Compute kernels using the canonical engine ----
    [kernels, kVars, kStats, hitStats, compStats] = computeSessionKernels(sessionData);

    % ---- Save noise matrices ----
    save(matrixFilePath, 'header', 'sideTypeNames', 'lr', ...
        'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', ...
        'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    fprintf('  Saved noise matrices: %s\n', matrixFilePath);

    % ---- Save kernels ----
    save(kernelFilePath, 'header', 'sideTypeNames', 'lr', ...
        'kernels', 'kVars', 'kStats', 'trialOutcomesAll', ...
        'changeSidesAll', 'changeIndicesAll', 'compStats', 'hitStats', '-v7.3');
    fprintf('  Saved kernels:        %s\n', kernelFilePath);

   % ---- Save kernel summary ----
    makeKernelSessionSummaries();

    % ---- Plot/export ----
    plotKernels(1, sprintf('%s (Probe %g%c)', baseName, probeDirDeg, char(176)), ...
    header, kernels, kVars, compStats, hitStats);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
    fprintf('  Saved plot:           %s\n', plotFilePath);

    % ---- Mark this probe direction as stale only after successful write ----
    staleProbeDirs(end+1) = probeDirDeg; %#ok<AGROW>
  end

  if numSkipped > 0
    fprintf(' makeKernels: Skipped %d previously processed files.\n', numSkipped);
  end
  staleProbeDirs = unique(staleProbeDirs);
end