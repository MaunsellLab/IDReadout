function looPredictor = computeIDQLOONoisePredictor(sessionAnalyses)
% computeIDQLOONoisePredictor
%
% Compute leave-one-session-out kernel-weighted signed-noise predictors.
% A common temporal kernel is estimated from the summed noise of all other
% sessions and then applied separately to each of the three direction streams.
%
% Outputs are aligned with the concatenated trial-table order used by
% makeIDQAcrossSessionSummary:
%   noisePredDir1, noisePredDir2, noisePredDir3

nSessions = numel(sessionAnalyses);

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
nStepFrames = numel(stepFrames);

predCell = cell(nSessions, 3);

sessionIndex = nan(nSessions, 1);
fileName = strings(nSessions, 1);
kernelStepSum = nan(nSessions, 1);
kernelStepMean = nan(nSessions, 1);
kernelStepMin = nan(nSessions, 1);
kernelStepMax = nan(nSessions, 1);
fracNegativeWeights = nan(nSessions, 1);
effectiveNFrames = nan(nSessions, 1);

weightsBySession = nan(nStepFrames, nSessions);

for iSession = 1:nSessions
    SA = sessionAnalyses{iSession};
    if numel(SA.tMS) ~= numel(tMS) || any(SA.tMS(:) ~= tMS(:))
        error('computeIDQLOONoisePredictor:TimeVectorMismatch', ...
            'Session %s has a different tMS vector.', SA.fileName);
    end
    if numel(SA.stepFrames) ~= numel(stepFrames) || any(SA.stepFrames(:) ~= stepFrames(:))
        error('computeIDQLOONoisePredictor:StepFrameMismatch', ...
            'Session %s has different stepFrames.', SA.fileName);
    end
    if size(SA.dirNoiseByFrameTrial, 3) ~= 3
        error('computeIDQLOONoisePredictor:ExpectedThreeDirs', ...
            'Session %s does not contain three direction streams.', SA.fileName);
    end

    looKernel = computeLOOKernel(sessionAnalyses, iSession);
    kStep = looKernel.meanDiff(stepFrames);
    kSum = sum(kStep, 'omitnan');

    if ~isfinite(kSum) || kSum == 0
        error('computeIDQLOONoisePredictor:BadKernelStepSum', ...
            'LOO kernel step sum is zero or nonfinite for session %s.', SA.fileName);
    end

    w = kStep(:) ./ kSum;

    for iDir = 1:3
        dirNoise = SA.dirNoiseByFrameTrial(stepFrames, :, iDir);
        predCell{iSession, iDir} = sum(dirNoise .* w, 1, 'omitnan')';
    end

    sessionIndex(iSession) = iSession;
    fileName(iSession) = string(SA.fileName);
    kernelStepSum(iSession) = kSum;
    kernelStepMean(iSession) = mean(kStep, 'omitnan');
    kernelStepMin(iSession) = min(kStep, [], 'omitnan');
    kernelStepMax(iSession) = max(kStep, [], 'omitnan');
    fracNegativeWeights(iSession) = mean(w < 0);
    effectiveNFrames(iSession) = 1 / sum(w .^ 2);
    weightsBySession(:, iSession) = w;
end

looPredictor = struct();
looPredictor.noisePredDir1 = vertcat(predCell{:, 1});
looPredictor.noisePredDir2 = vertcat(predCell{:, 2});
looPredictor.noisePredDir3 = vertcat(predCell{:, 3});
looPredictor.stepFrames = stepFrames;
looPredictor.tMS = tMS;
looPredictor.weightsBySession = weightsBySession;
looPredictor.weightDiagnostics = table( ...
    sessionIndex, fileName, kernelStepSum, kernelStepMean, kernelStepMin, ...
    kernelStepMax, fracNegativeWeights, effectiveNFrames);

end

%% ------------------------------------------------------------------------
function kernel = computeLOOKernel(sessionAnalyses, leaveOutSession)

allCorrectNoise = [];
allErrorNoise = [];

for iSession = 1:numel(sessionAnalyses)
    if iSession == leaveOutSession && numel(sessionAnalyses) > 1
        continue
    end
    SA = sessionAnalyses{iSession};
    T = SA.trialTable;

    idxUse = T.hasStepNoise;
    correctUse = idxUse & T.correct;
    errorUse = idxUse & ~T.correct;
    allCorrectNoise = [allCorrectNoise, SA.sumNoiseByFrameTrial(:, correctUse)]; %#ok<AGROW>
    allErrorNoise = [allErrorNoise, SA.sumNoiseByFrameTrial(:, errorUse)]; %#ok<AGROW>
end

kernel = struct();
kernel.meanCorrect = mean(allCorrectNoise, 2, 'omitnan');
kernel.meanError = mean(allErrorNoise, 2, 'omitnan');
kernel.meanDiff = kernel.meanCorrect - kernel.meanError;
kernel.nCorrect = size(allCorrectNoise, 2);
kernel.nError = size(allErrorNoise, 2);
kernel.nTrials = kernel.nCorrect + kernel.nError;

end
