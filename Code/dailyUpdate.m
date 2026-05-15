function dailyUpdate(replace, doBootstrap, nBoot, path, forceAcrossOffset)
% dailyUpdate  Run the daily MT-kernel analysis update pipeline.
%
%   dailyUpdate()
%       Runs with defaults:
%         replace     = false
%         doBootstrap = false
%         nBoot       = 100
%         path        = folderPath()
%
%   Pipeline:
%     1) convert any unconverted raw dat files
%     2) build missing/stale probe-specific kernel, noise-matrix, and plot outputs
%     3) build missing/stale probe-specific kernel-session summaries
%     4) refresh average-kernel plots for affected probe directions
%     5) refresh across-offset summaries/plots

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
  if nargin < 5 || isempty(forceAcrossOffset)
    forceAcrossOffset = false;
  end
  path = char(path);

  fprintf('>>> dailyUpdate start ---\n');

  % ---- Convert raw files ----
  fprintf('  >> convertIDRData start ---\n');
  convertIDRData;
  fprintf('  << convertIDRData complete ---\n');

  % ---- Session-level probe-specific kernels/noise matrices ----
  fprintf('  >> makeKernels start ---\n');
  [allProbeDirs, staleProbeDirs] = makeKernels(replace, path);
  anythingChanged = replace || ~isempty(staleProbeDirs);
  fprintf('  << makeKernels complete ---\n');

  allProbeDirs = unique(allProbeDirs);
  staleProbeDirs = unique(staleProbeDirs);

  % If replace is requested, refresh all known probe dirs.
  if replace
    refreshProbeDirs = allProbeDirs;
  else
    refreshProbeDirs = staleProbeDirs;
  end

  % ---- Per-session summaries ----
  if isempty(refreshProbeDirs)
    fprintf('    no stale probe-specific session outputs detected; skipping session summaries and averages.\n');
  else
    sessionsDirs = probeDirsToSessionDirs(refreshProbeDirs);

    fprintf('  >> makeKernelSessionSummaries start ---\n');
    makeKernelSessionSummaries('path', path, 'sessionsDirs', sessionsDirs, 'replace', replace, 'doBootstrap', false);
    fprintf('  << makeKernelSessionSummaries complete ---\n');

    % ---- Average-kernel plots ----
    fprintf('  >> kernelAverage start ---\n');
    for p = refreshProbeDirs(:).'
      probeTag = sprintf('probe%d', round(p));
      fprintf('      Updating average for %s\n', probeTag);
      kernelAverage(true, 50, 'dataFolder', fullfile(path, 'Data', probeTag, 'NoiseMatrices'), ...
          'plotFolder', fullfile(path, 'Plots', 'AverageKernels'), 'probeDirDeg', p);
    end
    fprintf('  << kernelAverage complete ---\n');
  end

  fprintf('  >> plotSideTypeKernelAverage start ---\n');
  plotSideTypeKernelAverage();
  fprintf('  << plotSideTypeKernelAverage complete ---\n');
  % ---- Across-offset summary update ----
  % This should run even when no single-session files were stale, because it
  % is cheap relative to the pipeline and keeps summary/plots synchronized
  % with any manual changes to summaries or exclusion rules.
  % ---- Across-offset summary update ----
  if anythingChanged || doBootstrap || forceAcrossOffset
    fprintf('  >> updateAcrossOffsetSummaries start ---\n');
    acrossOffsetSummary = updateAcrossOffsetSummaries(fullfile(path, 'Data'), ...
      'SaveFile', fullfile(path, 'Data', 'AcrossOffsetSummaries', 'IDR_acrossOffsetSummary.mat'), ...
      'PlotDir', fullfile(path, 'Plots', 'ReadoutFits'),'NBoot', nBoot, 'RandomSeed', 1, ...
      'MakePlots', true); %#ok<NASGU>
    fprintf('  << updateAcrossOffsetSummaries complete ---\n');
  else
    fprintf('    no session-level updates detected; skipping across-offset bootstrap/fits.\n');
  end
  fprintf('<<< Daily update complete ---\n');
end

%% probeDirsToSessionDirs()
function sessionsDirs = probeDirsToSessionDirs(probeDirs)

sessionsDirs = cell(1, numel(probeDirs));
for i = 1:numel(probeDirs)
  sessionsDirs{i} = sprintf('probe%d', round(probeDirs(i)));
end

end