function summaryTable = inspectIDQSessionAnalysis()
% inspectIDQSessionAnalysis
%
% Load IDQ session-analysis files and report compact diagnostics.
%
% Reads:
%   Data/SessionAnalysis/*_sessionAnalysis.mat
%
% Each file must contain:
%   sessionAnalysis
%
% Output:
%   summaryTable, also printed to command window.

domainPath = domainFolder(mfilename('fullpath'));

analysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
files = dir(fullfile(analysisFolder, '*_sessionAnalysis.mat'));

nFiles = numel(files);

fileName = strings(nFiles, 1);
nTrials = nan(nFiles, 1);
nStepNoiseTrials = nan(nFiles, 1);
noisyStepCoh = nan(nFiles, 1);

meanCorrect = nan(nFiles, 1);
meanCorrectStepNoise = nan(nFiles, 1);

meanRectNoise = nan(nFiles, 1);
stdRectNoise = nan(nFiles, 1);

kernelStepIntegral = nan(nFiles, 1);
kernelStepPeak = nan(nFiles, 1);
kernelStepMean = nan(nFiles, 1);

nCorrectKernel = nan(nFiles, 1);
nErrorKernel = nan(nFiles, 1);

for iFile = 1:nFiles

  filePath = fullfile(files(iFile).folder, files(iFile).name);
  load(filePath, 'sessionAnalysis');

  T = sessionAnalysis.trialTable;

  fileName(iFile) = string(sessionAnalysis.fileName);
  nTrials(iFile) = sessionAnalysis.nTrials;
  nStepNoiseTrials(iFile) = sessionAnalysis.nStepNoiseTrials;
  noisyStepCoh(iFile) = sessionAnalysis.noisyStepCoh;

  meanCorrect(iFile) = mean(T.correct, 'omitnan');

  stepNoiseIdx = T.hasStepNoise;
  meanCorrectStepNoise(iFile) = mean(T.correct(stepNoiseIdx), 'omitnan');

  x = sessionAnalysis.rectNoisePredictor(stepNoiseIdx);
  meanRectNoise(iFile) = mean(x, 'omitnan');
  stdRectNoise(iFile) = std(x, 'omitnan');

  k = sessionAnalysis.signedNoiseKernel.kernel;
  stepFrames = sessionAnalysis.stepFrames;

  kStep = k(stepFrames);

  kernelStepIntegral(iFile) = sum(kStep, 'omitnan');
  kernelStepMean(iFile) = mean(kStep, 'omitnan');
  kernelStepPeak(iFile) = max(abs(kStep), [], 'omitnan');

  nCorrectKernel(iFile) = sessionAnalysis.signedNoiseKernel.nCorrect;
  nErrorKernel(iFile) = sessionAnalysis.signedNoiseKernel.nError;

end

summaryTable = table( ...
  fileName, ...
  nTrials, ...
  nStepNoiseTrials, ...
  noisyStepCoh, ...
  meanCorrect, ...
  meanCorrectStepNoise, ...
  meanRectNoise, ...
  stdRectNoise, ...
  kernelStepIntegral, ...
  kernelStepMean, ...
  kernelStepPeak, ...
  nCorrectKernel, ...
  nErrorKernel);

disp(summaryTable);

end