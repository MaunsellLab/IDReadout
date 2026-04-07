function makeKernelSessionSummaries(varargin)
% makeKernelSessionSummaries
%
% Batch-create per-session tracking summaries from kernel files.

P = inputParser;
addParameter(P, 'kernelDir', '/Users/Shared/Data/IDReadout/Data/Kernels', @ischar);
addParameter(P, 'summaryDir', '/Users/Shared/Data/IDReadout/Data/KernelSummaries', @ischar);
addParameter(P, 'replace', false, @islogical);
addParameter(P, 'doBootstrap', true, @islogical);
addParameter(P, 'nBoot', 500, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(P, 'trackSideType', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(P, 'trackStepType', 2, @(x) isnumeric(x) && isscalar(x));
parse(P, varargin{:});
R = P.Results;

files = dir(fullfile(R.kernelDir, '*.mat'));

for i = 1:numel(files)
  kernelFile = fullfile(files(i).folder, files(i).name);

  try
    compileKernelSessionSummary(kernelFile, ...
      'summaryDir', R.summaryDir, ...
      'replace', R.replace, ...
      'doBootstrap', R.doBootstrap, ...
      'nBoot', R.nBoot, ...
      'trackSideType', R.trackSideType, ...
      'trackStepType', R.trackStepType);
  catch ME
    warning('makeKernelSessionSummaries:failed', ...
      'Failed on %s: %s', files(i).name, ME.message);
  end
end

end