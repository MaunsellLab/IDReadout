function directionDiagnostics = computeIDQDirectionDiagnostics(sessionAnalyses, trialTable)
% computeIDQDirectionDiagnostics
%
% Direction-specific diagnostics for IDQ across-session analysis.
%
% These are diagnostic only. The primary readout analysis remains collapsed
% across directions.

directionDiagnostics = struct();

directionDiagnostics.absolute = computeAbsoluteDirectionDiagnostics(sessionAnalyses, trialTable);
directionDiagnostics.aligned = computeAlignedDirectionDiagnostics(sessionAnalyses, trialTable);

end

%% ------------------------------------------------------------------------
function D = computeAbsoluteDirectionDiagnostics(sessionAnalyses, trialTable)

dirLabels = getAbsoluteDirectionLabels(sessionAnalyses);

D = struct();
D.mode = 'absolute';
D.modeLabel = 'Absolute directions';
D.directionLabels = dirLabels;

D.behavior = computeBehaviorByDriftDirection(trialTable, dirLabels);
D.kernels = computeAbsoluteNoiseKernels(sessionAnalyses, dirLabels);

end

%% ------------------------------------------------------------------------
function D = computeAlignedDirectionDiagnostics(sessionAnalyses, trialTable)

dirLabels = ["drift", "non-drift +1", "non-drift -1"];

D = struct();
D.mode = 'aligned';
D.modeLabel = 'Aligned to drift direction';
D.directionLabels = dirLabels;

% Behavior is not meaningful in aligned stream coordinates; every trial's
% selected drift direction is, by definition, the aligned "drift" direction.
D.behavior = table();

D.kernels = computeAlignedNoiseKernels(sessionAnalyses, dirLabels);

end

%% ------------------------------------------------------------------------
function dirLabels = getAbsoluteDirectionLabels(sessionAnalyses)

dirsDeg = sessionAnalyses{1}.sessionHeader.dirsDeg(:)';

for iSession = 2:numel(sessionAnalyses)
    thisDirsDeg = sessionAnalyses{iSession}.sessionHeader.dirsDeg(:)';
    if numel(thisDirsDeg) ~= numel(dirsDeg) || any(thisDirsDeg ~= dirsDeg)
        error('computeIDQDirectionDiagnostics:DirDegMismatch', 'sessionHeader.dirsDeg differs across sessions.');
    end
end

dirLabels = strings(1, numel(dirsDeg));
for iDir = 1:numel(dirsDeg)
    % dirLabels(iDir) = sprintf('%g deg', dirsDeg(iDir));
    dirLabels(iDir) = sprintf('%.0f°', dirsDeg(iDir));
end

end

%% ------------------------------------------------------------------------
function behavior = computeBehaviorByDriftDirection(trialTable, dirLabels)

nDirs = numel(dirLabels);

dirIndex = (1:nDirs)';
directionLabel = dirLabels(:);

nTrials = nan(nDirs, 1);
nCorrect = nan(nDirs, 1);
pCorrect = nan(nDirs, 1);
nStepNoiseTrials = nan(nDirs, 1);
pCorrectStepNoise = nan(nDirs, 1);
meanStepCoh = nan(nDirs, 1);
meanRectNoisePredictor = nan(nDirs, 1);

for iDir = 1:nDirs

    idx = trialTable.dirIndex == iDir;
    idxNoise = idx & trialTable.hasStepNoise;

    nTrials(iDir) = sum(idx);
    nCorrect(iDir) = sum(trialTable.correct(idx));
    pCorrect(iDir) = mean(trialTable.correct(idx), 'omitnan');

    nStepNoiseTrials(iDir) = sum(idxNoise);
    pCorrectStepNoise(iDir) = mean(trialTable.correct(idxNoise), 'omitnan');

    meanStepCoh(iDir) = mean(trialTable.stepCoh(idx), 'omitnan');

    if ismember('rectNoisePredictor', trialTable.Properties.VariableNames)
      if ismember('rectSumNoise', trialTable.Properties.VariableNames)
        meanRectNoisePredictor(iDir) = mean(trialTable.rectSumNoise(idxNoise), 'omitnan');
      elseif ismember('rectNoisePredictor', trialTable.Properties.VariableNames)
        meanRectNoisePredictor(iDir) = mean(trialTable.rectNoisePredictor(idxNoise), 'omitnan');
      end
    end
end

behavior = table( ...
    dirIndex, ...
    directionLabel, ...
    nTrials, ...
    nCorrect, ...
    pCorrect, ...
    nStepNoiseTrials, ...
    pCorrectStepNoise, ...
    meanStepCoh, ...
    meanRectNoisePredictor);

end

%% ------------------------------------------------------------------------
function kernels = computeAbsoluteNoiseKernels(sessionAnalyses, dirLabels)

nDirs = numel(dirLabels);
tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;

meanCorrect = nan(numel(tMS), nDirs);
meanError = nan(numel(tMS), nDirs);
meanDiff = nan(numel(tMS), nDirs);

nCorrect = nan(nDirs, 1);
nError = nan(nDirs, 1);
nTrials = nan(nDirs, 1);

for iDir = 1:nDirs

    allCorrect = [];
    allError = [];

    for iSession = 1:numel(sessionAnalyses)

        SA = sessionAnalyses{iSession};
        T = SA.trialTable;

        idxUse = T.hasStepNoise;

        sideIndex = T.sideIndex(:)';
        dirIndex = iDir;

        nTrial = height(T);
        noiseThisDir = nan(numel(tMS), nTrial);

        for iTrial = 1:nTrial
            noiseThisDir(:, iTrial) = squeeze(SA.noiseBySideDir( ...
                sideIndex(iTrial), dirIndex, :, iTrial));
        end

        allCorrect = [allCorrect, noiseThisDir(:, idxUse & T.correct)];
        allError = [allError, noiseThisDir(:, idxUse & ~T.correct)];
    end

    meanCorrect(:, iDir) = mean(allCorrect, 2, 'omitnan');
    meanError(:, iDir) = mean(allError, 2, 'omitnan');
    meanDiff(:, iDir) = meanCorrect(:, iDir) - meanError(:, iDir);

    nCorrect(iDir) = size(allCorrect, 2);
    nError(iDir) = size(allError, 2);
    nTrials(iDir) = nCorrect(iDir) + nError(iDir);
end

kernels = packageKernelTable(tMS, stepFrames, dirLabels, ...
    meanCorrect, meanError, meanDiff, nTrials, nCorrect, nError);

end

%% ------------------------------------------------------------------------
function kernels = computeAlignedNoiseKernels(sessionAnalyses, dirLabels)

nDirs = numel(dirLabels);
tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;

meanCorrect = nan(numel(tMS), nDirs);
meanError = nan(numel(tMS), nDirs);
meanDiff = nan(numel(tMS), nDirs);

nCorrect = nan(nDirs, 1);
nError = nan(nDirs, 1);
nTrials = nan(nDirs, 1);

for iRel = 1:nDirs

    allCorrect = [];
    allError = [];

    for iSession = 1:numel(sessionAnalyses)

        SA = sessionAnalyses{iSession};
        T = SA.trialTable;

        idxUse = T.hasStepNoise;

        sideIndex = T.sideIndex(:);
        driftDirIndex = T.dirIndex(:);

        nTrial = height(T);
        noiseThisRel = nan(numel(tMS), nTrial);

        for iTrial = 1:nTrial

            if iRel == 1
                dirIndex = driftDirIndex(iTrial);
            elseif iRel == 2
                dirIndex = mod(driftDirIndex(iTrial), 3) + 1;
            else
                dirIndex = mod(driftDirIndex(iTrial) - 2, 3) + 1;
            end

            noiseThisRel(:, iTrial) = squeeze(SA.noiseBySideDir( ...
                sideIndex(iTrial), dirIndex, :, iTrial));
        end

        allCorrect = [allCorrect, noiseThisRel(:, idxUse & T.correct)];
        allError = [allError, noiseThisRel(:, idxUse & ~T.correct)];
    end

    meanCorrect(:, iRel) = mean(allCorrect, 2, 'omitnan');
    meanError(:, iRel) = mean(allError, 2, 'omitnan');
    meanDiff(:, iRel) = meanCorrect(:, iRel) - meanError(:, iRel);

    nCorrect(iRel) = size(allCorrect, 2);
    nError(iRel) = size(allError, 2);
    nTrials(iRel) = nCorrect(iRel) + nError(iRel);
end

kernels = packageKernelTable(tMS, stepFrames, dirLabels, ...
    meanCorrect, meanError, meanDiff, nTrials, nCorrect, nError);

end

%% ------------------------------------------------------------------------
function kernels = packageKernelTable(tMS, stepFrames, dirLabels, ...
    meanCorrect, meanError, meanDiff, nTrials, nCorrect, nError)

nDirs = numel(dirLabels);

dirIndex = (1:nDirs)';
directionLabel = dirLabels(:);

stepMean = nan(nDirs, 1);
stepIntegral = nan(nDirs, 1);
stepPeakAbs = nan(nDirs, 1);

for iDir = 1:nDirs
    stepMean(iDir) = mean(meanDiff(stepFrames, iDir), 'omitnan');
    stepIntegral(iDir) = sum(meanDiff(stepFrames, iDir), 'omitnan');
    stepPeakAbs(iDir) = max(abs(meanDiff(stepFrames, iDir)), [], 'omitnan');
end

summaryTable = table( ...
    dirIndex, ...
    directionLabel, ...
    nTrials, ...
    nCorrect, ...
    nError, ...
    stepMean, ...
    stepIntegral, ...
    stepPeakAbs);

kernels = struct();
kernels.tMS = tMS;
kernels.stepFrames = stepFrames;
kernels.meanCorrect = meanCorrect;
kernels.meanError = meanError;
kernels.meanDiff = meanDiff;
kernels.summaryTable = summaryTable;

end

