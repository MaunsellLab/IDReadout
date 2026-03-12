function makeKernels(replace, path)
% NEW makeKernels  Generate kernel and plot files for all Lablib .mat data files.
%
%   makeKernels()
%       Uses folderPath() to get the base path.
%
%   makeKernels(path)
%       Uses the specified base path.
%
%   makeKernels(path, replace)
%       If replace is true, recomputes kernels and plots even if they exist.
%       If replace is false (default), only computes missing kernels/plots.
%
%   Directory conventions:
%       Data .mat files:      path/Data/*.mat   (excluding *_fileInfo.mat)
%       Kernel .mat files:    path/Kernels/<baseName>.mat
%       Kernel plot PDFs:     path/Plots/Kernels/<baseName>.pdf

  % ---- Handle inputs ----
  if nargin < 1 || isempty(replace)
    replace = false;
  end  
  if nargin < 2 || isempty(path)
    path = folderPath();   % user-defined helper
  end
  path = char(path);           % Ensure path is a char (fullfile is happier with chars on older MATLABs)

  % ---- Define directories and create as needed ----
  dataFolder   = fullfile(path, 'Data');
  kernelFolder = fullfile(path, 'Kernels');
  plotFolder   = fullfile(path, 'Plots', 'Kernels');
  if ~exist(kernelFolder, 'dir')
    mkdir(kernelFolder);
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
  names = {allMatFiles.name};                         % Exclude *_fileInfo.mat files
  isFileInfo = endsWith(names, '_fileInfo.mat');
  dataFiles = allMatFiles(~isFileInfo);
  if isempty(dataFiles)
    fprintf('No data .mat files (excluding *_fileInfo.mat) found in %s\n', dataFolder);
    return;
  end

  % ---- Process each data file ----
  for k = 1:numel(dataFiles)
    dataFileName = dataFiles(k).name;
    dataFilePath = fullfile(dataFolder, dataFileName);
    [~, baseName] = fileparts(dataFileName);
    kernelFilePath = fullfile(kernelFolder, [baseName '.mat']);
    plotFilePath   = fullfile(plotFolder,   [baseName '.pdf']);
    kernelExists = isfile(kernelFilePath);
    plotExists   = isfile(plotFilePath);
    if ~replace && kernelExists && plotExists
      fprintf('Skipping %s (kernel + plot already exist)\n', dataFileName);
      continue;
    end
    fprintf('Processing %s ...\n', dataFileName);

    % ---- Load header and trials from the data file ----
    S = load(dataFilePath, 'header', 'trials');
    header = S.header;
    trials = S.trials;

    % ---- Extract values ----
    frameRateHz = header.frameRateHz.data(1);
    msPerVFrame = 1000.0 / frameRateHz;

    % Sometimes header values are vectors; we only need the first entry
    preStepMS = header.preStepMS.data(1);
    stepMS    = header.stepMS.data(1);

    % Number of video frames up to (preStep + step)
    m = round((preStepMS + stepMS) / msPerVFrame);
    prefVals  = cellfun(@(t) t.trial.data.prefCohNoisePC,  trials);
    probeVals = cellfun(@(t) t.trial.data.probeCohNoisePC, trials);  
    prefCohNoisePC  = max(prefVals(prefVals  > 0));
    probeCohNoisePC = max(probeVals(probeVals > 0));

    % 2 x 2 x m kernel vectors (inc/dec, pref/probe, vFrame)
    kernels = nan(2, 2, m);
    % 2 x 2 kernel variance scalar (inc/dec, pref/probe)
    kVars = nan(2, 2);
    % 2 x 2 struct of sufficient statistics (inc/dec, pref/probe)
    kStats = repmat(struct('nCorrect',0,'nWrong',0,'sumCorrect',[],'sumWrong',[],'sigma2',nan), 2, 2);    
    % 1 x 2 cell of trial outcome vectors (inc/dec); pref & probe always same
    trialOutcomes = cell(1, 2);
    nHits   = nan(1, 2);
    nTrials = nan(1, 2);

    % ---- Compute kernels ----
    % sideType == 0 means difference evidence (change - noChange)
    sideType = 0;
    if sideType == 0
      ampScale = sqrt(2);
    else
      ampScale = 1;
    end
    for s = 1:2
      [prefMat, probeMat, trialOutcomes{s}] = extractNoiseMatrices(header, trials, s, sideType);
      save (sprintf('Matrices%d', s), 'prefMat', 'probeMat', 'trialOutcomes');
    
      nTrials(s) = numel(trialOutcomes{s});
      nHits(s)   = nTrials(s) - sum(trialOutcomes{s});
    
      [kernels(s,1,:), kVars(s,1), kStats(s,1)] = meanPsychKernel(prefMat,  trialOutcomes{s}, ampScale*prefCohNoisePC);
      [kernels(s,2,:), kVars(s,2), kStats(s,2)] = meanPsychKernel(probeMat, trialOutcomes{s}, ampScale*probeCohNoisePC);
    end
    [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame);

    % ---- Plot kernels ----
    % Use figure 1 (as before);
    plotKernels(1, baseName, header, kernels, kVars, kIntegrals, R, RVar, nHits, nTrials);
    if ~exist(plotFolder, 'dir')  % Ensure plot directory exists (in case deleted it mid-run)
      mkdir(plotFolder);
    end
    print(gcf, plotFilePath, '-dpdf');

    % ---- Save kernel data ----
    if ~exist(kernelFolder, 'dir')
      mkdir(kernelFolder);
    end
    save(kernelFilePath, 'header', 'kernels', 'kVars', 'kStats', 'trialOutcomes', ...
                         'kIntegrals', 'R', 'RVar', 'nHits', 'nTrials', '-v7.3');
    fprintf('  Saved kernels: %s\n', kernelFilePath);
    fprintf('  Saved plot:    %s\n', plotFilePath);
  end
end