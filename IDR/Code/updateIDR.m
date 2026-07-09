function updateIDR
% dailyUpdate  Run the daily MT-kernel analysis update pipeline.
%%
%   Pipeline:
%     1) convert any unconverted raw dat files
%     2) build missing/stale probe-specific kernel, noise-matrix, and plot outputs
%     3) build missing/stale probe-specific kernel-session summaries
%     4) refresh average-kernel plots for affected probe directions
%     5) refresh across-offset summaries/plots

fprintf('>>> dailyUpdate start\n');
% replace = false;
nBoot = 500;
% animal = 'All';
% animal = 'Meetz';
% animal = 'Neesha';

% ---- Convert raw files ----
fprintf('  >> convertIDRData start\n');
convertIDRData;
fprintf('  << convertIDRData complete\n');

% ---- Session-level probe-specific data files ----
fprintf('  >> makeProbeSessions start\n');
makeProbeSessions(); 
fprintf('  << makeProbeSessions complete\n');

% ---- Make probe session kernels and kernel plots ----
fprintf('  >> makeKernels start\n');
makeKernels(); 
fprintf('  << makeKernels complete\n');

% ---- Make beta session summaries ----
fprintf('  >> makeBetaSessionData start\n');
makeBetaSessionData();
fprintf('  << makeBetaSessionData complete\n');

% The following stages of the update all involve across-session analyses of
% various sorts (including session beta regressions, which use an
% across-session kernel for temporally weighting the noise signal.


% ---- Probe-session regressions and across-offset beta summary ----

%IDQ Analysis

fprintf('  >> makeBetaKernel start\n');
makeBetaKernel('Animal', 'Meetz');
fprintf('  >> makeBetaKernel complete\n');

fprintf('  >> makeIDRGainAcrossSessionPrefFit start\n');
makeIDRGainAcrossSessionPrefFit('Animal', 'Meetz');
fprintf('  >> makeIDRGainAcrossSessionPrefFit complete\n');

fprintf('  >> makeIDRGainProbeOffsetFits start\n');
makeIDRGainProbeOffsetFits('Animal', 'Meetz');
fprintf('  >> makeIDRGainProbeOffsetFits complete\n');

fprintf('  >> makeIDRGainProbeOffsetSummary start\n');
makeIDRGainProbeOffsetSummary('Animal', 'Meetz');
fprintf('  >> makeIDRGainProbeOffsetSummary complete\n');

%IDR Analysis

% ---- Average Kernels ----
fprintf('  >> kernelAverage start\n');
for p = [10, 25, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 100, 'probeDirDeg', p, 'Verbose', true, 'Animal', 'Neesha');
end
fprintf('  << kernelAverage complete\n');

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 25, 45, 90, 135, 179], 'Animal', 'Neesha');
fprintf('  << plotSideTypeKernelAverage complete\n');

% ---- Across-offset summary update ----
fprintf('  >> updateAcrossOffsetSummaries start\n');
updateAcrossOffsetSummaries([], 'NBoot', nBoot, 'RandomSeed', 1, 'Animal', 'Neesha');
fprintf('  << updateAcrossOffsetSummaries complete\n');

fprintf('<<< daily update complete\n');
end
