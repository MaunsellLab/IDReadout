function sessionAnalysisAll = makeIDQSessionAnalyses(processedFile)
% makeIDQSessionAnalyses
%
% Build compact session-level analysis products from IDQ processed session
% files. These outputs are intended to support both daily session PDFs and
% across-session analyses.
%
% Usage:
%   makeIDQSessionAnalyses()
%   makeIDQSessionAnalyses(processedFile)
%
% Input:
%   processedFile: optional full path to one processed session .mat file.
%                  If omitted, all .mat files in Data/ProcessedSessions are
%                  processed.
%
% Output files:
%   Data/SessionAnalysis/<sessionName>_sessionAnalysis.mat
%
% Required variables in processed session files:
%   header
%   sessionHeader
%   trialData
%   noiseBySideDir
%
% Required trialData fields:
%   correct
%   validIdx
%   sideIndex
%   chosenSideIndex
%   dirIndex
%   stepCoh
%   hasStepNoise
%
% Required sessionHeader fields:
%   fileName
%   frameRateHz
%   preStepMS
%   stepMS

if nargin < 1
  processedFile = '';
end

thisFolder = domainFolder(mfilename('fullpath'));

processedFolder = fullfile(thisFolder, 'Data', 'ProcessedSessions');
analysisFolder  = fullfile(thisFolder, 'Data', 'SessionAnalysis');

if ~exist(analysisFolder, 'dir')
  mkdir(analysisFolder);
end

if isempty(processedFile)
  files = dir(fullfile(processedFolder, '*.mat'));
  processedFiles = fullfile({files.folder}, {files.name})';
else
  processedFiles = {processedFile};
end

sessionAnalysisAll = cell(numel(processedFiles), 1);

for iFile = 1:numel(processedFiles)
  sessionAnalysis = makeOneSessionAnalysis(processedFiles{iFile}, analysisFolder);
  sessionAnalysisAll{iFile} = sessionAnalysis;
end

if isscalar(sessionAnalysisAll)
  sessionAnalysisAll = sessionAnalysisAll{1};
end

end

%% ------------------------------------------------------------------------
function sessionAnalysis = makeOneSessionAnalysis(processedFile, analysisFolder)

load(processedFile, 'header', 'sessionHeader', 'trialData', 'noiseBySideDir');

sessionName = sessionHeader.fileName;
sessionName = erase(sessionName, '.mat');
sessionName = erase(sessionName, '.dat');

fprintf('makeIDQSessionAnalyses: %s\n', sessionName);

[nSides, nDirs, nFrames, nTrialsNoise] = size(noiseBySideDir);

nTrials = numel(trialData.correct);

if nSides ~= 2
  error('Expected noiseBySideDir side dimension to be 2.');
end

if nDirs ~= 3
  error('Expected noiseBySideDir direction dimension to be 3.');
end

if nTrialsNoise ~= nTrials
  error('noiseBySideDir trial dimension does not match trialData.');
end

frameMS = 1000 / sessionHeader.frameRateHz;
tMS = (0:nFrames - 1)' * frameMS;
stepFrames = find(tMS >= sessionHeader.preStepMS & ...
                  tMS <  sessionHeader.preStepMS + sessionHeader.stepMS);
noiseMeasures = computeNoiseMeasuresByFrameTrial(noiseBySideDir, trialData);

sumNoiseByFrameTrial = noiseMeasures.sumNoiseByFrameTrial;
meanNoiseByFrameTrial = noiseMeasures.meanNoiseByFrameTrial;
dirNoiseByFrameTrial = noiseMeasures.dirNoiseByFrameTrial;

% Primary flat-readout predictor: changed-side dir1 + dir2 + dir3, averaged over the step window.
rectSumNoise = mean(sumNoiseByFrameTrial(stepFrames, :), 1, 'omitnan')';

% Convenience/display predictor.  If used in a gain fit, flat readout predicts gain = 3.
rectMeanNoise = mean(meanNoiseByFrameTrial(stepFrames, :), 1, 'omitnan')';

% Legacy diagnostic only: old drift-minus-nondrift convention.
% rectDriftMinusNonNoise = mean(driftMinusNonNoiseByFrameTrial(stepFrames, :), 1, 'omitnan')';

% Compatibility alias for downstream code.  This is now the summed-noise predictor.
rectNoisePredictor = rectSumNoise;

sumNoiseKernel = computeNoiseKernel( ...
    sumNoiseByFrameTrial, trialData.correct, trialData.hasStepNoise);

meanNoiseKernel = computeNoiseKernel( ...
    meanNoiseByFrameTrial, trialData.correct, trialData.hasStepNoise);

% driftMinusNonNoiseKernel = computeNoiseKernel( ...
    % driftMinusNonNoiseByFrameTrial, trialData.correct, trialData.hasStepNoise);
noisyStepCoh = unique(trialData.stepCoh(trialData.hasStepNoise));
if numel(noisyStepCoh) ~= 1
  error('Expected exactly one noisy step coherence in session %s.', sessionName);
end

sessionAnalysis = struct();

trialTable = makeTrialTable(trialData, sessionName);
trialTable.rectNoisePredictor = rectNoisePredictor;                 % compatibility alias
trialTable.rectSumNoise = rectSumNoise;                             % primary
trialTable.rectMeanNoise = rectMeanNoise;                           % display / gain=3 convention
% trialTable.rectDriftMinusNonNoise = rectDriftMinusNonNoise;         % legacy diagnostic
sessionAnalysis.trialTable = trialTable;

sessionAnalysis.fileName = sessionName;
sessionAnalysis.sourceProcessedFile = processedFile;
sessionAnalysis.createdBy = mfilename;
sessionAnalysis.createdAt = datetime('now');

sessionAnalysis.header = header;
sessionAnalysis.sessionHeader = sessionHeader;

sessionAnalysis.nTrials = nTrials;
sessionAnalysis.nStepNoiseTrials = sum(trialData.hasStepNoise);
sessionAnalysis.noisyStepCoh = noisyStepCoh;

sessionAnalysis.tMS = tMS;
sessionAnalysis.stepFrames = stepFrames;

sessionAnalysis.trialTable = trialTable;
sessionAnalysis.primaryNoisePredictor = 'rectSumNoise';
sessionAnalysis.primaryNoiseDefinition = ...
  'changed-side dir1 + dir2 + dir3 noise, averaged over step frames';

% Primary and diagnostic framewise measures.
sessionAnalysis.sumNoiseByFrameTrial = sumNoiseByFrameTrial;
sessionAnalysis.meanNoiseByFrameTrial = meanNoiseByFrameTrial;
% sessionAnalysis.driftMinusNonNoiseByFrameTrial = driftMinusNonNoiseByFrameTrial;
sessionAnalysis.dirNoiseByFrameTrial = dirNoiseByFrameTrial;

% Primary and diagnostic rectangular predictors.
sessionAnalysis.rectNoisePredictor = rectNoisePredictor;   % compatibility alias
sessionAnalysis.rectSumNoise = rectSumNoise;
sessionAnalysis.rectMeanNoise = rectMeanNoise;
% sessionAnalysis.rectDriftMinusNonNoise = rectDriftMinusNonNoise;

% Primary and diagnostic kernels.
sessionAnalysis.sumNoiseKernel = sumNoiseKernel;
sessionAnalysis.meanNoiseKernel = meanNoiseKernel;
% sessionAnalysis.driftMinusNonNoiseKernel = driftMinusNonNoiseKernel;

% Compatibility alias for downstream code.  This is now summed-noise, not
% drift-minus-nondrift signed noise.
sessionAnalysis.signedNoiseByFrameTrial = sumNoiseByFrameTrial;
sessionAnalysis.signedNoiseKernel = sumNoiseKernel;
sessionAnalysis.noiseBySideDir = noiseBySideDir;

outFile = fullfile(analysisFolder, sprintf('%s_sessionAnalysis.mat', sessionName));
save(outFile, 'sessionAnalysis', '-v7.3');

% fprintf('  saved %s\n', outFile);

end

%% ------------------------------------------------------------------------
function trialTable = makeTrialTable(trialData, sessionName)

nTrials = numel(trialData.correct);

trialTable = table();

trialTable.sessionName = repmat(string(sessionName), nTrials, 1);
trialTable.trialIndex = (1:nTrials)';

trialTable.correct = logical(trialData.correct(:));
trialTable.validIdx = trialData.validIdx(:);

trialTable.sideIndex = trialData.sideIndex(:);
trialTable.chosenSideIndex = trialData.chosenSideIndex(:);
trialTable.dirIndex = trialData.dirIndex(:);

trialTable.stepCoh = trialData.stepCoh(:);
trialTable.hasStepNoise = logical(trialData.hasStepNoise(:));

end

%% ------------------------------------------------------------------------
function noiseMeasures = computeNoiseMeasuresByFrameTrial(noiseBySideDir, trialData)

[~, nDirs, nFrames, nTrials] = size(noiseBySideDir);

if nDirs ~= 3
  error('computeNoiseMeasuresByFrameTrial:ExpectedThreeDirs', ...
    'Expected noiseBySideDir to have 3 direction streams.');
end

sumNoiseByFrameTrial = nan(nFrames, nTrials);
meanNoiseByFrameTrial = nan(nFrames, nTrials);
dirNoiseByFrameTrial = nan(nFrames, nTrials, nDirs);

for iTrial = 1:nTrials

  sideIndex = trialData.sideIndex(iTrial);
  driftDir = trialData.dirIndex(iTrial);

  if ~isfinite(sideIndex) || ~isfinite(driftDir)
    continue
  end

  sideIndex = double(sideIndex);
  driftDir = double(driftDir);

  if sideIndex < 1 || sideIndex > size(noiseBySideDir, 1) || ...
      driftDir < 1 || driftDir > nDirs
    continue
  end

  % Changed-side noise for all three direction streams.
  %
  % noiseBySideDir is:
  %   side x direction x frame x trial
  %
  % changedSideDirNoise is:
  %   direction x frame
  changedSideDirNoise = squeeze(noiseBySideDir(sideIndex, :, :, iTrial));

  if size(changedSideDirNoise, 1) ~= nDirs
    changedSideDirNoise = changedSideDirNoise';
  end

  dirNoiseByFrameTrial(:, iTrial, :) = changedSideDirNoise';

  % Primary flat-readout measure:
  % all three changed-side direction streams contribute with the same sign.
  sumNoiseByFrameTrial(:, iTrial) = sum(changedSideDirNoise, 1)';
  meanNoiseByFrameTrial(:, iTrial) = mean(changedSideDirNoise, 1)';

end

noiseMeasures = struct();
noiseMeasures.sumNoiseByFrameTrial = sumNoiseByFrameTrial;
noiseMeasures.meanNoiseByFrameTrial = meanNoiseByFrameTrial;
noiseMeasures.dirNoiseByFrameTrial = dirNoiseByFrameTrial;

end

%% ------------------------------------------------------------------------
function kernel = computeNoiseKernel(noiseByFrameTrial, correct, hasStepNoise)

correct = logical(correct(:));
hasStepNoise = logical(hasStepNoise(:));

useTrials = hasStepNoise & isfinite(correct);

correctTrials = useTrials & correct;
errorTrials = useTrials & ~correct;

meanCorrect = mean(noiseByFrameTrial(:, correctTrials), 2, 'omitnan');
meanError = mean(noiseByFrameTrial(:, errorTrials), 2, 'omitnan');
kernel = struct();

kernel.kernel = meanCorrect - meanError;
kernel.meanCorrect = meanCorrect;
kernel.meanError = meanError;

kernel.nTrials = sum(useTrials);
kernel.nCorrect = sum(correctTrials);
kernel.nError = sum(errorTrials);

end