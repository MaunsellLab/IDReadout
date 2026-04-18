function dailyUpdate(replace, doBootstrap, nBoot, path)
% dailyUpdate  Run the daily MT-kernel analysis update pipeline.
%
%   dailyUpdate()
%       Runs with defaults:
%         replace     = false
%         doBootstrap = false
%         nBoot       = 100
%         path        = folderPath()
%
%   dailyUpdate(replace, doBootstrap, nBoot, path)
%       replace:
%         If true, recompute all eligible session outputs.
%         If false, only compute missing session outputs.
%
%       doBootstrap:
%         If true, run hierarchical bootstrap in kernelAverageByProbe.
%
%       nBoot:
%         Number of bootstrap repetitions.
%
%       path:
%         Base project folder.
%
%   Pipeline:
%     1) build any missing/stale session kernel, matrix, and plot outputs
%     2) collect probe directions affected by those updates
%     3) refresh per-probe averages only for affected probe directions
%     4) refresh across-offset summaries/plots

  % ---- Handle inputs ----
  if nargin < 1 || isempty(replace)
    replace = false;
  end
  if nargin < 2 || isempty(doBootstrap)
    doBootstrap = false;
  end
  if nargin < 3 || isempty(nBoot)
    nBoot = 100;
  end
  if nargin < 4 || isempty(path)
    path = folderPath();
  end
  path = char(path);

  fprintf('--- Daily update start ---\n');
  convertIDRData;

  % ---- Session-level update ----
  staleProbeDirs = makeKernels(replace, path);

  if isempty(staleProbeDirs) && ~doBootstrap
    fprintf(' No session-level updates detected. Skipping per-probe averages.\n');
  else
    if doBootstrap
      staleProbeDirs = [45, 180];
    end
    fprintf('Refreshing probe averages for: %s\n', mat2str(staleProbeDirs));
    for p = staleProbeDirs(:).'
      fprintf(' Updating average for probe %g\n', p);
      kernelAverageByProbe(path, p, doBootstrap, nBoot);
    end
  end

  % ---- Across-offset summary update ----
  % Keep this call unconditional if you want cross-offset summaries to stay
  % synchronized with any existing per-probe summary files.
  fprintf(' Updating across-offset summaries\n');
  acrossOffsetSummary = updateAcrossOffsetSummaries('/Users/Shared/Data/IDReadout/Data/KernelSummaries', ...
        'NBoot', 10, 'MakePlots', false, 'ScaleSideType', 'change', 'ScaleStepType', 'inc'); 

  fprintf('--- Daily update complete ---\n');

end