
%% makeProbeSessions
clc;
% cleanupObj = initProjectPath(); %#ok<NASGU>
fprintf('  >> makeProbeSessions start\n');
makeProbeSessions(false);
fprintf('  << makeProbeSessions complete\n');

%% makeKernels

% Make the kernels
fprintf('  >> makeKernels start\n');
makeKernels(false);
fprintf('  << makeKernels complete\n');

%% ---- Average-kernel plots ----
clc
fprintf('  >> kernelAverage start\n');
for p = [10, 25, 45, 90, 135, 179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 1, 'probeDirDeg', p, 'FileSelectionArgs', {'Bin179With180', true});
end
fprintf('  << kernelAverage complete\n');

%% ---- Average-kernel plots ----
clc
fprintf('  >> kernelAverage start\n');
for p = [179, 180]
  fprintf('      updating average for probe %d\n', p);
  kernelAverage(true, 5, 'probeDirDeg', p, 'FileSelectionArgs', {'Bin179With180', true});
end
fprintf('  << kernelAverage complete\n');

%% plotSideTypeKernelAverage
% Summary plot of kernels from increment/change side

% cleanupObj = initProjectPath();
fprintf('  >> plotSideTypeKernelAverage start\n');
plotSideTypeKernelAverage('ProbeDirs', [10, 25, 45, 90, 135, 179]);
fprintf('  << plotSideTypeKernelAverage complete\n');

%% updateAcrossOffsetSummaries
clc;
acrossOffsetSummary = updateAcrossOffsetSummaries([], 'Verbose', true, 'NBoot', 1000, 'RandomSeed', 1, ...
    'FileSelectionArgs', {'Bin179With180', true});

%% updateScalarNoiseRegression
clc;
batch = updateScalarNoiseRegression(false);

%% updateScalarNoiseRegression
clc;
summaryTable = summarizeScalarNoiseRegression(false);

