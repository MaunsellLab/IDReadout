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

if nargin < 1
  processedFile = '';
end
thisFolder = domainFolder(mfilename('fullpath'));

processedFolder = fullfile(thisFolder, 'Data', 'ProcessedSessions');
analysisFolder  = validFolder(fullfile(thisFolder, 'Data', 'SessionAnalysis'));

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

[~, ~, nFrames, nTrialsNoise] = size(noiseBySideDir);
nTrials = numel(trialData.correct);
if nTrialsNoise ~= nTrials
  error('noiseBySideDir trial dimension does not match trialData.');
end

frameMS = 1000 / sessionHeader.frameRateHz;
tMS = (0:nFrames - 1)' * frameMS;
stepFrames = find(tMS >= sessionHeader.preStepMS & tMS <  sessionHeader.preStepMS + sessionHeader.stepMS);
noiseMeasures = computeNoiseMeasuresByFrameTrial(noiseBySideDir, trialData);
changeSumNoiseByFrameTrial = noiseMeasures.changeSumNoiseByFrameTrial;
changeDirNoiseByFrameTrial = noiseMeasures.changeDirNoiseByFrameTrial;
noChangeSumNoiseByFrameTrial = noiseMeasures.noChangeSumNoiseByFrameTrial;
noChangeDirNoiseByFrameTrial = noiseMeasures.noChangeDirNoiseByFrameTrial;

% Rectangular predictors can be computed within each session. Retain explicit
% names for both temporal windows; plotIDQSummary copies the selected family
% into the neutral noisePredDir1-3 variables expected by the gain fitter.
fullFrames = (1:nFrames)';
rectStepChangeNoisePredDir1 = mean(changeDirNoiseByFrameTrial(stepFrames, :, 1), 1, 'omitnan')';
rectStepChangeNoisePredDir2 = mean(changeDirNoiseByFrameTrial(stepFrames, :, 2), 1, 'omitnan')';
rectStepChangeNoisePredDir3 = mean(changeDirNoiseByFrameTrial(stepFrames, :, 3), 1, 'omitnan')';
rectFullChangeNoisePredDir1 = mean(changeDirNoiseByFrameTrial(fullFrames, :, 1), 1, 'omitnan')';
rectFullChangeNoisePredDir2 = mean(changeDirNoiseByFrameTrial(fullFrames, :, 2), 1, 'omitnan')';
rectFullChangeNoisePredDir3 = mean(changeDirNoiseByFrameTrial(fullFrames, :, 3), 1, 'omitnan')';

rectStepNoChangeNoisePredDir1 = mean(noChangeDirNoiseByFrameTrial(stepFrames, :, 1), 1, 'omitnan')';
rectStepNoChangeNoisePredDir2 = mean(noChangeDirNoiseByFrameTrial(stepFrames, :, 2), 1, 'omitnan')';
rectStepNoChangeNoisePredDir3 = mean(noChangeDirNoiseByFrameTrial(stepFrames, :, 3), 1, 'omitnan')';
rectFullNoChangeNoisePredDir1 = mean(noChangeDirNoiseByFrameTrial(fullFrames, :, 1), 1, 'omitnan')';
rectFullNoChangeNoisePredDir2 = mean(noChangeDirNoiseByFrameTrial(fullFrames, :, 2), 1, 'omitnan')';
rectFullNoChangeNoisePredDir3 = mean(noChangeDirNoiseByFrameTrial(fullFrames, :, 3), 1, 'omitnan')';

changeSumNoiseKernel = computeNoiseKernel(changeSumNoiseByFrameTrial, ...
  trialData.correct, trialData.hasStepNoise);
noChangeSumNoiseKernel = computeNoiseKernel(noChangeSumNoiseByFrameTrial, ...
  trialData.correct, trialData.hasStepNoise);
noisyStepCoh = unique(trialData.stepCoh(trialData.hasStepNoise));
if numel(noisyStepCoh) ~= 1
  error('Expected exactly one noisy step coherence in session %s.', sessionName);
end

sessionAnalysis = struct();

trialTable = makeTrialTable(trialData, sessionName);
trialTable.rectStepChangeNoisePredDir1 = rectStepChangeNoisePredDir1;
trialTable.rectStepChangeNoisePredDir2 = rectStepChangeNoisePredDir2;
trialTable.rectStepChangeNoisePredDir3 = rectStepChangeNoisePredDir3;
trialTable.rectFullChangeNoisePredDir1 = rectFullChangeNoisePredDir1;
trialTable.rectFullChangeNoisePredDir2 = rectFullChangeNoisePredDir2;
trialTable.rectFullChangeNoisePredDir3 = rectFullChangeNoisePredDir3;

trialTable.rectStepNoChangeNoisePredDir1 = rectStepNoChangeNoisePredDir1;
trialTable.rectStepNoChangeNoisePredDir2 = rectStepNoChangeNoisePredDir2;
trialTable.rectStepNoChangeNoisePredDir3 = rectStepNoChangeNoisePredDir3;
trialTable.rectFullNoChangeNoisePredDir1 = rectFullNoChangeNoisePredDir1;
trialTable.rectFullNoChangeNoisePredDir2 = rectFullNoChangeNoisePredDir2;
trialTable.rectFullNoChangeNoisePredDir3 = rectFullNoChangeNoisePredDir3;

% Backward-compatible aliases: the original unqualified predictors are the
% changed-side predictors.
trialTable.rectStepNoisePredDir1 = rectStepChangeNoisePredDir1;
trialTable.rectStepNoisePredDir2 = rectStepChangeNoisePredDir2;
trialTable.rectStepNoisePredDir3 = rectStepChangeNoisePredDir3;
trialTable.rectFullNoisePredDir1 = rectFullChangeNoisePredDir1;
trialTable.rectFullNoisePredDir2 = rectFullChangeNoisePredDir2;
trialTable.rectFullNoisePredDir3 = rectFullChangeNoisePredDir3;
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
sessionAnalysis.rectangularPredictorDefinitions = { ...
  'change step: changed-side direction noise averaged over step frames'; ...
  'change full: changed-side direction noise averaged over all trial frames'; ...
  'no-change step: no-change-side direction noise averaged over step frames'; ...
  'no-change full: no-change-side direction noise averaged over all trial frames'; ...
  'all predictor signs retain the sign of physical directional coherence'};
sessionAnalysis.relativeDirectionDefinition = ...
  ['For both sides, drift-relative direction is defined from trialTable.dirIndex; ' ...
   'the no-change-side drift-relative stream is the absolute direction matching ' ...
   'the change-side drift direction.'];

% Primary and diagnostic frame-wise measures.
sessionAnalysis.changeSumNoiseByFrameTrial = changeSumNoiseByFrameTrial;
sessionAnalysis.changeDirNoiseByFrameTrial = changeDirNoiseByFrameTrial;
sessionAnalysis.noChangeSumNoiseByFrameTrial = noChangeSumNoiseByFrameTrial;
sessionAnalysis.noChangeDirNoiseByFrameTrial = noChangeDirNoiseByFrameTrial;

% Backward-compatible aliases: the original unqualified measures are from
% the changed side.
sessionAnalysis.sumNoiseByFrameTrial = changeSumNoiseByFrameTrial;
sessionAnalysis.dirNoiseByFrameTrial = changeDirNoiseByFrameTrial;

% Primary and diagnostic kernels.
sessionAnalysis.changeSumNoiseKernel = changeSumNoiseKernel;
sessionAnalysis.noChangeSumNoiseKernel = noChangeSumNoiseKernel;
sessionAnalysis.sumNoiseKernel = changeSumNoiseKernel;

sessionAnalysis.noiseBySideDir = noiseBySideDir;

outFile = fullfile(analysisFolder, sprintf('%s_sessionAnalysis.mat', sessionName));
save(outFile, 'sessionAnalysis', '-v7.3');

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
trialTable.changeSideIndex = trialTable.sideIndex;
trialTable.noChangeSideIndex = 3 - trialTable.sideIndex;
trialTable.chosenSideIndex = trialData.chosenSideIndex(:);
trialTable.dirIndex = trialData.dirIndex(:);


trialTable.stepCoh = trialData.stepCoh(:);
trialTable.hasStepNoise = logical(trialData.hasStepNoise(:));

end

%% ------------------------------------------------------------------------
function noiseMeasures = computeNoiseMeasuresByFrameTrial(noiseBySideDir, trialData)

[nSides, nDirs, nFrames, nTrials] = size(noiseBySideDir);

if nDirs ~= 3
  error('computeNoiseMeasuresByFrameTrial:ExpectedThreeDirs', ...
    'Expected noiseBySideDir to have 3 direction streams.');
end
if nSides ~= 2
  error('computeNoiseMeasuresByFrameTrial:ExpectedTwoSides', ...
    'Expected noiseBySideDir to have 2 sides.');
end

changeSumNoiseByFrameTrial = nan(nFrames, nTrials);
changeDirNoiseByFrameTrial = nan(nFrames, nTrials, nDirs);
noChangeSumNoiseByFrameTrial = nan(nFrames, nTrials);
noChangeDirNoiseByFrameTrial = nan(nFrames, nTrials, nDirs);

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

  noChangeSideIndex = nSides + 1 - sideIndex;

  % Noise for all three direction streams on the changed and no-change sides.
  %
  % noiseBySideDir is:
  %   side x direction x frame x trial
  %
  % Each extracted side matrix is:
  %   direction x frame
  changeSideDirNoise = squeeze(noiseBySideDir(sideIndex, :, :, iTrial));
  noChangeSideDirNoise = squeeze(noiseBySideDir(noChangeSideIndex, :, :, iTrial));

  if size(changeSideDirNoise, 1) ~= nDirs
    changeSideDirNoise = changeSideDirNoise';
  end
  if size(noChangeSideDirNoise, 1) ~= nDirs
    noChangeSideDirNoise = noChangeSideDirNoise';
  end

  changeDirNoiseByFrameTrial(:, iTrial, :) = changeSideDirNoise';
  noChangeDirNoiseByFrameTrial(:, iTrial, :) = noChangeSideDirNoise';

  changeSumNoiseByFrameTrial(:, iTrial) = sum(changeSideDirNoise, 1)';
  noChangeSumNoiseByFrameTrial(:, iTrial) = sum(noChangeSideDirNoise, 1)';
end

noiseMeasures = struct();
noiseMeasures.changeSumNoiseByFrameTrial = changeSumNoiseByFrameTrial;
noiseMeasures.changeDirNoiseByFrameTrial = changeDirNoiseByFrameTrial;
noiseMeasures.noChangeSumNoiseByFrameTrial = noChangeSumNoiseByFrameTrial;
noiseMeasures.noChangeDirNoiseByFrameTrial = noChangeDirNoiseByFrameTrial;

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
