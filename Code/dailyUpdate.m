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
  fprintf('---  convertIDRData start ---\n');
  convertIDRData;
  fprintf('---  convertIDRData complete ---\n');

  % ---- Session-level update ----
  fprintf('---  makeKernels start ---\n');
  [allProbeDirs, staleProbeDirs] = makeKernels(replace, path);
  fprintf('---  makeKernels complete ---\n');

  % ---- Save kernel summary ----
  fprintf('---  makeKernelSessionSummaries start ---\n');
  makeKernelSessionSummaries();
  fprintf('---  makeKernelSessionSummaries complete ---\n');

  if isempty(staleProbeDirs) && ~doBootstrap
    fprintf('--- No session-level updates detected; skipping per-probe averages.\n');
  else
    if doBootstrap
      staleProbeDirs = allProbeDirs;
    end
    for p = staleProbeDirs(:).'
      fprintf('--- Updating average for probe %g\n', p);
      fprintf('---  kernelAverageByProbe start ---\n');
      kernelAverageByProbe(path, p, doBootstrap, nBoot);
      fprintf('---  kernelAverageByProbe complete ---\n');
    end
  end

  % ---- Across-offset summary update ----
  fprintf('---  updateAcrossOffsetSummaries start ---\n');
  acrossOffsetSummary = updateAcrossOffsetSummaries([], 'NBoot', nBoot, 'RandomSeed', 1); %#ok<NASGU>
  fprintf('---  updateAcrossOffsetSummaries complete ---\n');
  fprintf('--- Daily update complete ---\n');
end