
%% makeKernels

% Make the kernels
fprintf('  >> makeKernels start\n');
makeKernels(true);
fprintf('  << makeKernels complete\n');

%% ---- Average-kernel plots for just the selected probes ----

fprintf('  >> kernelAverage start\n');
p = 179;
fprintf('      updating average for probe%d\n', p);
kernelAverage(true, 5, 'probeDirDeg', p, 'Verbose', true, 'FileSelectionArgs', {'Bin179With180', true});
fprintf('  << kernelAverage complete\n');
%% ---- Average-kernel plots ----

fprintf('  >> kernelAverage start\n');
for p = [10, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe%d\n', p);
  kernelAverage(true, 5, 'probeDirDeg', p, 'FileSelectionArgs', {'Bin179With180', true});
end
fprintf('  << kernelAverage complete\n');

%% plotSideTypeKernelAverage
% Summary plot of kernels from increment/change side

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 45, 90, 135, 179]);
fprintf('  << plotSideTypeKernelAverage complete\n');

%% updateAcrossOffsetSummaries

acrossOffsetSummary = updateAcrossOffsetSummaries([], 'Verbose', true, 'NBoot', 5, 'RandomSeed', 1, ...
    'FileSelectionArgs', {'Bin179With180', true});

