
%% makeProbeSessions
clc;
fprintf('  >> makeProbeSessions start\n');
makeProbeSessions(true);
fprintf('  << makeProbeSessions complete\n');

%% makeKernels

% Make the kernels
fprintf('  >> makeKernels start\n');
makeKernels(true);
fprintf('  << makeKernels complete\n');

%% ---- Average-kernel plots ----
clc
fprintf('  >> kernelAverage start\n');
for p = [10, 25, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 500, 'probeDirDeg', p, 'SummarySideType', 'Chosen', 'FileSelectionArgs', {'Bin179With180', true});
end
fprintf('  << kernelAverage complete\n');

%% ---- Average-kernel plots ----
clc
fprintf('  >> kernelAverage start\n');
for p = [179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 5, 'probeDirDeg', p, 'SummarySideType', 'Chosen', 'FileSelectionArgs', {'Bin179With180', true});
end
fprintf('  << kernelAverage complete\n');

%% plotSideTypeKernelAverage
% Summary plot of kernels from increment/change side

fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 25, 45, 90, 135, 179], 'SideType', 'Chosen');
fprintf('  << plotSideTypeKernelAverage complete\n');

%% updateAcrossOffsetSummaries
clc;
acrossOffsetSummary = updateAcrossOffsetSummaries([], 'Verbose', true, 'NBoot', 5, 'RandomSeed', 1, ...
    'FileSelectionArgs', {'Bin179With180', true});

