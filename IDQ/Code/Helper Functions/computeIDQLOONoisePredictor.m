function looPredictor = computeIDQLOONoisePredictor(sessionAnalyses)
% computeIDQLOONoisePredictor
%
% Compute leave-one-session-out kernel-weighted signed-noise predictors.
%
% Compute leave-one-session-out kernel-weighted summed-noise predictors.
%
% For each session:
%   1. Build the across-session summed-noise kernel excluding that session.
%   2. Use its step-period values as weights.
%   3. Normalize weights by their sum.
%   4. Apply weights to that session's sumNoiseByFrameTrial.
%
% Output:
%   looPredictor.xNoiseLOO is one vector aligned with the concatenated
%   trialTable order used by makeIDQAcrossSessionSummary.

nSessions = numel(sessionAnalyses);

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
nStepFrames = numel(stepFrames);

xNoiseLOOCell = cell(nSessions, 1);

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

    looKernel = computeLOOKernel(sessionAnalyses, iSession);
    kStep = looKernel.meanDiff(stepFrames);

    kSum = sum(kStep, 'omitnan');

    if ~isfinite(kSum) || kSum == 0
        error('computeIDQLOONoisePredictor:BadKernelStepSum', ...
            'LOO kernel step sum is zero or nonfinite for session %s.', SA.fileName);
    end

    w = kStep(:) ./ kSum;

    xNoiseLOO = sum(SA.sumNoiseByFrameTrial(stepFrames, :) .* w, 1, 'omitnan')';
    xNoiseLOOCell{iSession} = xNoiseLOO;

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

looPredictor.xNoiseLOO = vertcat(xNoiseLOOCell{:});
looPredictor.stepFrames = stepFrames;
looPredictor.tMS = tMS;
looPredictor.weightsBySession = weightsBySession;

looPredictor.weightDiagnostics = table( ...
    sessionIndex, ...
    fileName, ...
    kernelStepSum, ...
    kernelStepMean, ...
    kernelStepMin, ...
    kernelStepMax, ...
    fracNegativeWeights, ...
    effectiveNFrames);

end

%% ------------------------------------------------------------------------
function kernel = computeLOOKernel(sessionAnalyses, leaveOutSession)
% computeLOOKernel: computer a leave-one-out kernel.  Returns a single 
% kernels if only one session is passed in.

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
  allCorrectNoise = [allCorrectNoise, SA.sumNoiseByFrameTrial(:, correctUse)];
  allErrorNoise = [allErrorNoise, SA.sumNoiseByFrameTrial(:, errorUse)];
end

meanCorrect = mean(allCorrectNoise, 2, 'omitnan');
meanError = mean(allErrorNoise, 2, 'omitnan');

kernel = struct();
kernel.meanCorrect = meanCorrect;
kernel.meanError = meanError;
kernel.meanDiff = meanCorrect - meanError;
kernel.nCorrect = size(allCorrectNoise, 2);
kernel.nError = size(allErrorNoise, 2);
kernel.nTrials = kernel.nCorrect + kernel.nError;

end