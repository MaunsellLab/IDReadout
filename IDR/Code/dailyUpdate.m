function dailyUpdate()
% dailyUpdate  Run the daily MT-kernel analysis update pipeline.
%%
%   Pipeline:
%     1) convert any unconverted raw dat files
%     2) build missing/stale probe-specific kernel, noise-matrix, and plot outputs
%     3) build missing/stale probe-specific kernel-session summaries
%     4) refresh average-kernel plots for affected probe directions
%     5) refresh across-offset summaries/plots

fprintf('>>> dailyUpdate start\n');
replace = false;
doKernelBootstrap = false;
% doBetaBootstrap = false;
nKernelBoot = 500;
% nBetaBoot = 500;
% bootstrapSeed = 1;

% ---- Convert raw files ----
fprintf('  >> convertIDRData start\n');
convertIDRData;
fprintf('  << convertIDRData complete\n');

% ---- Session-level probe-specific data files ----
fprintf('  >> makeProbeSessions start\n');
[allProbeDirs, staleProbeDirs] = makeProbeSessions(replace); %#ok<ASGLU>
fprintf('  << makeProbeSessions complete\n');

% ---- Beta temporal weights ----
fprintf('  >> beta weight update start\n');
[nBetaCreated, ~] = makeBetaSessionData(false);
weightsPath = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'FullSessions', 'BetaAnalysis', ...
  'AcrossSessions', 'BetaWeights.mat');
betaWeightsChanged = nBetaCreated > 0 || ~isfile(weightsPath);
if betaWeightsChanged
  [~, issueTable] = validateBetaSessionData();
  if ~isempty(issueTable)
    error('dailyUpdate:BetaValidationFailed', 'validateBetaSessionData reported %d issue(s).', height(issueTable));
  end
  makeBetaKernel();
  makeBetaWeights();
end

plotAcrossOffsetBetaSummary('MakeRatioPlot', true, 'MakeReadoutPlot', true, 'MakeDiagnosticPlots', true);
fprintf('  << beta weight update complete\n');

% ---- Make probe session kernels and kernel plots ----
fprintf('  >> makeKernels start\n');
[allKernelDirs, staleKernelDirs] = makeKernels(replace); %#ok<ASGLU>
anythingChanged = replace || ~isempty(staleProbeDirs) || ~isempty(staleKernelDirs);
if replace
  refreshProbeDirs = unique([allProbeDirs, allKernelDirs]); %#ok<UNRCH>
else
  refreshProbeDirs = unique([staleProbeDirs, staleKernelDirs]);
end
fprintf('  << makeKernels complete\n');

% ---- Average Kernels ----
if isempty(refreshProbeDirs)
  fprintf('       no stale probe-specific session outputs detected; skipping session summaries and averages.\n');
else
  fprintf('  >> kernelAverage start\n');
  for p = refreshProbeDirs(:).'
    fprintf('      updating average for probe %d\n', p);
    kernelAverage(true, 100, 'probeDirDeg', p);
  end
  fprintf('  << kernelAverage complete\n');
end

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 25, 45, 90, 135, 179]);
fprintf('  << plotSideTypeKernelAverage complete\n');

% ---- Across-offset summary update ----
% This should run even when no single-session files were stale, because it is cheap relative to the pipeline and 
% keeps summary/plots synchronized with manual changes to summaries or exclusion rules.
if anythingChanged || doKernelBootstrap
  fprintf('  >> updateAcrossOffsetSummaries start\n');
  acrossOffsetSummary = updateAcrossOffsetSummaries([], 'NBoot', nKernelBoot, 'RandomSeed', 1, ...
    'FileSelectionArgs', {'Bin179With180', true}); %#ok<NASGU>
  fprintf('  << updateAcrossOffsetSummaries complete\n');
else
  fprintf('      no session-level updates detected; skipping across-offset bootstrap/fits.\n');
end

% ---- Kernel-versus-beta comparison ----
kernelSummaryPath = fullfile(domainFolder(mfilename('fullpath')), 'Data', 'AcrossOffsetSummaries', ...
  'IDR_acrossOffsetSummary.mat');
betaSummaryPath = fullfile(domainFolder(mfilename('fullpath')),  'Data', 'AcrossOffsetSummaries', ...
  'IDR_acrossOffsetBetaSummary.mat');
if isfile(kernelSummaryPath) && isfile(betaSummaryPath) &&  (anythingChanged || betaWeightsChanged || doKernelBootstrap)
  fprintf('  >> kernel-beta comparison plot start\n');
  plotKernelBetaReadoutComparison();
  fprintf('  << kernel-beta comparison plot complete\n');
end

fprintf('<<< daily update complete\n');
end