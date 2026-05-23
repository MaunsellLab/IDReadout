
%% makeKernels

% Make the kernels
fprintf('  >> makeKernels start\n');
makeKernels(false);
fprintf('  << makeKernels complete\n');

%% ---- Average-kernel plots ----

fprintf('  >> kernelAverage start\n');
% for p = [1, 10, 45, 90, 135, 179, 180]
for p = 45
    fprintf('      updating average for probe%d\n', p);
  kernelAverage(true, 5, 'probeDirDeg', p, ...
    'FileSelectionArgs', {'MaxParentNProbeDirections', 1});
end
fprintf('  << kernelAverage complete\n');%% ---- Average-kernel plots ----
%% ---- Average-kernel plots ----

fprintf('  >> kernelAverage start\n');
for p = [1, 10, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe%d\n', p);
  kernelAverage(true, 5, 'probeDirDeg', p);
end
fprintf('  << kernelAverage complete\n');

%% plotSideTypeKernelAverage
% Summary plot of kernels from increment/change side

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage();
fprintf('  << plotSideTypeKernelAverage complete\n');

%% updateAcrossOffsetSummaries

acrossOffsetSummary = updateAcrossOffsetSummaries([], 'Verbose', true, 'NBoot', 5, 'RandomSeed', 1);

