function makeStrategyModels(replace, path)
% makeStrategyModels  Generate cached behavioral strategy model files
% for all Lablib .mat data files.
%
%   makeStrategyModels()
%       Uses folderPath() to get the base path.
%
%   makeStrategyModels(replace)
%       If replace is true, recomputes model files even if they exist.
%       If replace is false (default), only computes missing model files.
%
%   makeStrategyModels(replace, path)
%       Uses the specified base path.
%
%   Directory conventions:
%       Data .mat files:      path/Data/*.mat   (excluding *_fileInfo.mat)
%       Model .mat files:     path/Models/<baseName>.mat
%
%   Each model file contains:
%       M   struct returned by strategyForensics(dataFilePath)

  % ---- Handle inputs ----
  if nargin < 1 || isempty(replace)
    replace = false;
  end
  if nargin < 2 || isempty(path)
    path = folderPath();   % user-defined helper
  end
  if isstring(path)        % Ensure path is char for compatibility
    path = char(path);
  end

  % ---- Define directories and create as needed ----
  dataFolder  = fullfile(path, 'Data');
  modelFolder = fullfile(path, 'Models');
  if ~exist(modelFolder, 'dir')
    mkdir(modelFolder);
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
  nDone = 0;
  nSkip = 0;
  nFail = 0;

  for k = 1:numel(dataFiles)
    dataFileName = dataFiles(k).name;
    dataFilePath = fullfile(dataFolder, dataFileName);
    [~, baseName] = fileparts(dataFileName);

    modelFilePath = fullfile(modelFolder, [baseName '.mat']);
    modelExists = isfile(modelFilePath);

    if ~replace && modelExists
      fprintf('Skipping %s (model already exists)\n', dataFileName);
      nSkip = nSkip + 1;
      continue;
    end

    fprintf('Processing %s ...\n', dataFileName);

    try
      % ---- Run per-session behavioral strategy analysis ----
      M = strategyForensics(dataFilePath);

      % ---- Save cached model file ----
      if ~exist(modelFolder, 'dir')
        mkdir(modelFolder);
      end
      save(modelFilePath, 'M', '-v7.3');

      fprintf('  Saved model: %s\n', modelFilePath);
      nDone = nDone + 1;

    catch ME
      fprintf(2, '  FAILED: %s\n', dataFileName);
      fprintf(2, '  %s\n', ME.message);
      nFail = nFail + 1;
    end
  end

  % ---- Summary ----
  fprintf('\nmakeStrategyModels complete\n');
  fprintf('  Computed: %d\n', nDone);
  fprintf('  Skipped:  %d\n', nSkip);
  fprintf('  Failed:   %d\n', nFail);
end