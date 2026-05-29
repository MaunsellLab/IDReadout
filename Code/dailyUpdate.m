function dailyUpdate()
% dailyUpdate  Run the daily MT-kernel analysis update pipeline.
%%
%   Pipeline:
%     1) convert any unconverted raw dat files
%     2) build missing/stale probe-specific kernel, noise-matrix, and plot outputs
%     3) build missing/stale probe-specific kernel-session summaries
%     4) refresh average-kernel plots for affected probe directions
%     5) refresh across-offset summaries/plots

  replace = false;
  doBootstrap = false;
  nBoot = 500;

  fprintf('>>> dailyUpdate start\n');

  % ---- Convert raw files ----
  fprintf('  >> convertIDRData start\n');
  convertIDRData;
  fprintf('  << convertIDRData complete\n');

  % ---- Session-level probe-specific kernels/noise matrices ----
  fprintf('  >> makeKernels start\n');
  [allProbeDirs, staleProbeDirs] = makeKernels(replace);
  anythingChanged = replace || ~isempty(staleProbeDirs);
  fprintf('  << makeKernels complete\n');

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
    fprintf('      no stale probe-specific session outputs detected; skipping session summaries and averages.\n');
  else
    sessionsDirs = probeDirsToSessionDirs(refreshProbeDirs);

    fprintf('  >> makeKernelSessionSummaries start\n');
    makeKernelSessionSummaries('sessionsDirs', sessionsDirs, 'replace', replace);
    fprintf('  << makeKernelSessionSummaries complete\n');

    % ---- Average-kernel plots ----
    fprintf('  >> kernelAverage start\n');
    for p = refreshProbeDirs(:).'
      probeTag = sprintf('probe%d', round(p));
      fprintf('      updating average for %s\n', probeTag);
      kernelAverage(true, 50, 'probeDirDeg', p);
    end
    fprintf('  << kernelAverage complete\n');
  end

  fprintf('  >> plotSideTypeKernelAverage start\n');
  plotSideTypeKernelAverage();
  fprintf('  << plotSideTypeKernelAverage complete\n');

  % ---- Across-offset summary update ----
  % This should run even when no single-session files were stale, because it
  % is cheap relative to the pipeline and keeps summary/plots synchronized
  % with any manual changes to summaries or exclusion rules.
  % ---- Across-offset summary update ----
  if anythingChanged || doBootstrap
    fprintf('  >> updateAcrossOffsetSummaries start\n');
    acrossOffsetSummary = updateAcrossOffsetSummaries([], 'Verbose', true, 'NBoot', 5, 'RandomSeed', 1, ...
      'FileSelectionArgs', {'Bin179With180', true}); %#ok<NASGU>
    % acrossOffsetSummary = updateAcrossOffsetSummaries([], 'NBoot', 5, 'RandomSeed', 1); %#ok<NASGU>
    fprintf('  << updateAcrossOffsetSummaries complete\n');
  else
    fprintf('      no session-level updates detected; skipping across-offset bootstrap/fits.\n');
  end
  fprintf('<<< daily update complete\n');
end

%% probeDirsToSessionDirs()
function sessionsDirs = probeDirsToSessionDirs(probeDirs)

sessionsDirs = cell(1, numel(probeDirs));
for i = 1:numel(probeDirs)
  sessionsDirs{i} = sprintf('probe%d', round(probeDirs(i)));
end

end