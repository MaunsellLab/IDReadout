function makeKernelSessionSummaries(varargin)
% makeKernelSessionSummaries
%
% % Batch utility for rebuilding KernelSummary files from existing kernel files.
% Routine makeKernels runs create KernelSummary files automatically.
% Use this function for backfilling, refreshing summaries after format changes,
% or computing optional per-session bootstrap diagnostics.

P = inputParser;
addParameter(P, 'kernelDir', '/Users/Shared/Data/IDReadout/Data/Kernels', @ischar);
addParameter(P, 'summaryDir', '/Users/Shared/Data/IDReadout/Data/KernelSummaries', @ischar);
addParameter(P, 'sessionsDirs', {}, @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(P, 'replace', false, @islogical);
addParameter(P, 'doBootstrap', false, @islogical);
addParameter(P, 'nBoot', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(P, 'trackSideType', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'trackStepType', 2, @(x) isnumeric(x) && isscalar(x));
parse(P, varargin{:});
R = P.Results;
R.path = folderPath();

sessionsDirs = normalizeSessionsDirs(R.sessionsDirs, R.path);

if isempty(sessionsDirs)
  % ---- Backward-compatible flat-folder mode ----
  processOneKernelDir(R.kernelDir, R.summaryDir, R);
else
  % ---- New probe-folder mode ----
  for s = 1:numel(sessionsDirs)
    sessionDir = sessionsDirs{s};

    kernelDir = fullfile(sessionDir, 'Kernels');
    summaryDir = fullfile(sessionDir, 'KernelSummaries');

    if ~exist(kernelDir, 'dir')
      warning('makeKernelSessionSummaries:MissingKernelDir', ...
        'Kernel directory not found: %s', kernelDir);
      continue;
    end

    if ~exist(summaryDir, 'dir')
      mkdir(summaryDir);
    end

    processOneKernelDir(kernelDir, summaryDir, R);
  end
end
end

%% processOneKernelDir()
function processOneKernelDir(kernelDir, summaryDir, R)

if ~exist(summaryDir, 'dir')
  mkdir(summaryDir);
end

files = dir(fullfile(kernelDir, '*.mat'));
for i = 1:numel(files)
  kernelFile = fullfile(files(i).folder, files(i).name);
  try
    compileKernelSessionSummary(kernelFile, ...
      'summaryDir', summaryDir, ...
      'replace', R.replace, ...
      'doBootstrap', R.doBootstrap, ...
      'nBoot', R.nBoot, ...
      'trackSideType', R.trackSideType, ...
      'trackStepType', R.trackStepType);
  catch ME
    warning('makeKernelSessionSummaries:FailedFile', ...
      'Failed on %s: %s', files(i).name, ME.message);
  end
end
end

%% normalizeSessionsDirs()
function sessionsDirs = normalizeSessionsDirs(sessionsDirsIn, pathRoot)

if isempty(sessionsDirsIn)
  sessionsDirs = {};
  return;
end

if ischar(sessionsDirsIn) || isstring(sessionsDirsIn)
  sessionsDirsIn = cellstr(sessionsDirsIn);
end

sessionsDirs = cell(size(sessionsDirsIn));

for i = 1:numel(sessionsDirsIn)
  d = char(sessionsDirsIn{i});

  % If user passed 'probe45', interpret relative to Data/.
  if ~isfolder(d)
    candidate = fullfile(pathRoot, 'Data', d);
    if isfolder(candidate)
      d = candidate;
    end
  end

  sessionsDirs{i} = d;
end
end