function makeKernels(replace, path)
% makeKernels  Generate kernel and plot files for all converted session files.
%
%   makeKernels()
%       Uses folderPath() to get the base path.
%
%   makeKernels(replace)
%       If replace is true, recomputes kernels and plots even if they exist.
%       If replace is false (default), only computes missing kernels/plots.
%
%   makeKernels(replace, path)
%       Uses the specified base path.
%
%   Directory conventions:
%       Converted data files:   path/Data/Converted/*.mat   (excluding *_fileInfo.mat)
%       Kernel .mat files:      path/Data/Kernels/<baseName>.mat
%       Noise matrix .mat:      path/Data/NoiseMatrices/<baseName>.mat
%       Kernel plot PDFs:       path/Plots/Kernels/<baseName>.pdf

  % ---- Handle inputs ----
  if nargin < 1 || isempty(replace)
    replace = false;
  end
  if nargin < 2 || isempty(path)
    path = folderPath();
  end
  path = char(path);   % older MATLABs are happier with char for fullfile

  % ---- Define directories and create as needed ----
  dataFolder   = fullfile(path, 'Data', 'Converted');
  kernelFolder = fullfile(path, 'Data', 'Kernels');
  matrixFolder = fullfile(path, 'Data', 'NoiseMatrices');
  plotFolder   = fullfile(path, 'Plots', 'Kernels');

  if ~exist(kernelFolder, 'dir')
    mkdir(kernelFolder);
  end
  if ~exist(matrixFolder, 'dir')
    mkdir(matrixFolder);
  end
  if ~exist(plotFolder, 'dir')
    mkdir(plotFolder);
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

  sideTypeNames = {'diff', 'change', 'noChange', 'L', 'R', 'RF', 'Opp'};

  % ---- Process each data file ----
  for k = 1:numel(dataFiles)
    dataFileName = dataFiles(k).name;
    dataFilePath = fullfile(dataFolder, dataFileName);
    [~, baseName] = fileparts(dataFileName);

    kernelFilePath = fullfile(kernelFolder, [baseName '.mat']);
    matrixFilePath = fullfile(matrixFolder, [baseName '.mat']);
    plotFilePath   = fullfile(plotFolder,   [baseName '.pdf']);

    kernelExists = isfile(kernelFilePath);
    matrixExists = isfile(matrixFilePath);
    plotExists   = isfile(plotFilePath);

    if ~replace && kernelExists && matrixExists && plotExists
      fprintf('Skipping %s (kernel + matrix + plot already exist)\n', dataFileName);
      continue;
    end

    fprintf('Processing %s ...\n', dataFileName);

    % ---- Load header and trials from the converted data file ----
    S = load(dataFilePath, 'header', 'trials');
    header = S.header;
    trials = S.trials;

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
    fprintf('  Saved kernels: %s\n', kernelFilePath);

    % ---- Plot/export ----
    plotKernels(1, baseName, header, kernels, kVars, compStats, hitStats);
    exportgraphics(gcf, plotFilePath, 'ContentType', 'vector');
    fprintf('  Saved plot:    %s\n', plotFilePath);
  end
end