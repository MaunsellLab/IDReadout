function updateIDR()
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
doKernelBootstrap = true;
nBoot = 500;
animal = 'Neesha';

% ---- Convert raw files ----
fprintf('  >> convertIDRData start\n');
convertIDRData;
fprintf('  << convertIDRData complete\n');

% ---- Session-level probe-specific data files ----
fprintf('  >> makeProbeSessions start\n');
staleProbeDirs = makeProbeSessions('Animal', animal, 'Replace', false); 
fprintf('  << makeProbeSessions complete\n');

% ---- Make probe session kernels and kernel plots ----
fprintf('  >> makeKernels start\n');
makeKernels('Animal', animal, 'Replace', false); 
fprintf('  << makeKernels complete\n');

% ---- Make beta session summaries ----
fprintf('  >> makeBetaSessionData start\n');
makeBetaSessionData('Animal', animal, 'Replace', false);
fprintf('  << makeBetaSessionData complete\n');

% The following stages of the update all involve across-session analyses of
% various sorts (including session beta regressions, which use an
% across-session kernel for temporally weighting the noise signal.

% ---- Beta temporal weights ----
fprintf('  >> beta fits start\n');
makeBetaKernel('Animal', animal);
fprintf('  << beta fits complete\n');

% ---- Probe-session regressions and across-offset beta summary ----
fprintf('  >> make, fit and plot regressions start\n');
makeRegressions('Replace', replace, 'Animal', animal);
fitAcrossOffsetBetaMeasurements('NBoot', nBoot, 'Animal', animal);
plotAcrossOffsetBetaSummary('Animal', animal);
fprintf('  << make, fit and plot regressions complete\n');

% ---- Average Kernels ----
fprintf('  >> kernelAverage start\n');
for p = [10, 25, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 100, 'probeDirDeg', p, 'Verbose', true, 'Animal', animal);
end
fprintf('  << kernelAverage complete\n');

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 25, 45, 90, 135, 179], 'Animal', animal);
fprintf('  << plotSideTypeKernelAverage complete\n');

% ---- Across-offset summary update ----
% This should run even when no single-session files were stale, because it is cheap relative to the pipeline and 
% keeps summary/plots synchronized with manual changes to summaries or exclusion rules.
if ~isempty(staleProbeDirs) || doKernelBootstrap
  fprintf('  >> updateAcrossOffsetSummaries start\n');
  acrossOffsetSummary = updateAcrossOffsetSummaries([], 'NBoot', nBoot, 'RandomSeed', 1, 'Animal', animal, ...
    'Bin179With180', true); %#ok<NASGU>
  fprintf('  << updateAcrossOffsetSummaries complete\n');
else
  fprintf('      no session-level updates detected; skipping across-offset bootstrap/fits.\n');
end

fprintf('<<< daily update complete\n');
end