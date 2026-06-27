function weightData = makeBetaWeights()
% makeBetaWeights  Create normalized all-session and leave-one-out weights.
%
% Reads:
%   Data/AcrossOffsetSummaries/BetaKernel_AllSessions.mat
%
% Saves:
%   Data/AcrossOffsetSummaries/BetaWeights.mat
%
% The weighting function is the unsmoothed preferred-noise kernel over the
% 250-ms step interval. Weights are normalized to sum to 1, so a constant
% +1% coherence-noise waveform produces +1% effective coherence.
%
% Leave-one-session-out kernels are reconstructed exactly from the saved
% per-session correct/error means and trial counts.

baseFolder = domainFolder(mfilename('fullpath'));
acrossFolder = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries');
kernelPath = fullfile(acrossFolder, 'BetaKernel_AllSessions.mat');
S = load(kernelPath, 'kernelData');
K = S.kernelData;
requiredFields = { 'tMS', 'stepOnsetMS', 'stepMS', 'sessionFileNames', 'sessionMeanCorrectPC', ...
  'sessionMeanErrorPC', 'sessionNCorrect', 'sessionNError'};
missing = requiredFields(~isfield(K, requiredFields));
if ~isempty(missing)
  error('makeBetaWeights:MissingFields', ...
    'kernelData is missing: %s', strjoin(missing, ', '));
end
tMS = double(K.tMS(:)');
stepStartMS = double(K.stepOnsetMS);
stepEndMS = stepStartMS + double(K.stepMS);
stepMask = tMS >= stepStartMS & tMS < stepEndMS;
stepTMS = tMS(stepMask);
nSessions = numel(K.sessionFileNames);
meanCorrect = K.sessionMeanCorrectPC;
meanError = K.sessionMeanErrorPC;
nCorrect = double(K.sessionNCorrect(:));
nError = double(K.sessionNError(:));

if size(meanCorrect,1) ~= nSessions || size(meanError,1) ~= nSessions
  error('makeBetaWeights:SessionCountMismatch', 'Per-session mean arrays do not match sessionFileNames.');
end
if numel(nCorrect) ~= nSessions || numel(nError) ~= nSessions
  error('makeBetaWeights:CountMismatch', 'Per-session trial counts do not match sessionFileNames.');
end
if any(nCorrect <= 0) || any(nError <= 0)
  error('makeBetaWeights:EmptyOutcomeClass', 'Every session must contain both correct and error trials.');
end

% ---- All-session pooled kernel ----
allKernelPC = pooledKernel(meanCorrect, meanError, nCorrect, nError, true(nSessions,1));
allStepKernelPC = allKernelPC(stepMask);
allKernelSum = sum(allStepKernelPC);

assertUsableKernelSum(allKernelSum, 'all-session');

allWeights = allStepKernelPC / allKernelSum;

% ---- Leave-one-session-out kernels and weights ----
looKernelPC = nan(nSessions, numel(tMS));
looStepKernelPC = nan(nSessions, numel(stepTMS));
looWeights = nan(nSessions, numel(stepTMS));
looKernelSum = nan(nSessions, 1);
looNCorrect = nan(nSessions, 1);
looNError = nan(nSessions, 1);

for iSession = 1:nSessions
  includeMask = true(nSessions,1);
  includeMask(iSession) = false;

  thisKernel = pooledKernel(meanCorrect, meanError, ...
    nCorrect, nError, includeMask);
  thisStepKernel = thisKernel(stepMask);
  thisSum = sum(thisStepKernel);

  assertUsableKernelSum(thisSum, sprintf('leave-one-out session %d (%s)', iSession, K.sessionFileNames{iSession}));

  looKernelPC(iSession,:) = thisKernel;
  looStepKernelPC(iSession,:) = thisStepKernel;
  looWeights(iSession,:) = thisStepKernel / thisSum;
  looKernelSum(iSession) = thisSum;
  looNCorrect(iSession) = sum(nCorrect(includeMask));
  looNError(iSession) = sum(nError(includeMask));
end

% ---- Output structure ----
weightData = struct();
weightData.version = 1;
weightData.method = ...
  ['unsmoothed preferred-noise kernel over the step interval, ' ...
   'normalized so weights sum to 1'];

weightData.timeOrigin = K.timeOrigin;
weightData.stepOnsetMS = stepStartMS;
weightData.stepMS = K.stepMS;
weightData.tMS = tMS;
weightData.stepMask = stepMask;
weightData.stepTMS = stepTMS;

weightData.allSessionKernelPC = allKernelPC;
weightData.allSessionStepKernelPC = allStepKernelPC;
weightData.allSessionKernelSum = allKernelSum;
weightData.allSessionWeights = allWeights;
weightData.allSessionWeightSum = sum(allWeights);

weightData.sessionFileNames = K.sessionFileNames;
weightData.leaveOneOutKernelPC = looKernelPC;
weightData.leaveOneOutStepKernelPC = looStepKernelPC;
weightData.leaveOneOutKernelSum = looKernelSum;
weightData.leaveOneOutWeights = looWeights;
weightData.leaveOneOutWeightSums = sum(looWeights, 2);
weightData.leaveOneOutNCorrect = looNCorrect;
weightData.leaveOneOutNError = looNError;

weightData.nSessions = nSessions;
weightData.createdBy = mfilename;
weightData.createdDate = datetime('now');

% ---- Final consistency checks ----
if abs(weightData.allSessionWeightSum - 1) > 1e-12
  error('makeBetaWeights:AllWeightNormalizationFailed', 'All-session weights do not sum to 1.');
end

if any(abs(weightData.leaveOneOutWeightSums - 1) > 1e-12)
  error('makeBetaWeights:LOOWeightNormalizationFailed', 'At least one leave-one-out weight vector does not sum to 1.');
end

outputPath = fullfile(acrossFolder, 'BetaWeights.mat');
save(outputPath, 'weightData', '-v7.3');
end

% -------------------------------------------------------------------------
function kernelPC = pooledKernel(meanCorrect, meanError, ...
  nCorrect, nError, includeMask)

totalCorrect = sum(nCorrect(includeMask));
totalError = sum(nError(includeMask));
if totalCorrect <= 0 || totalError <= 0
  error('makeBetaWeights:NoOutcomeTrials', 'Pooled data must contain both correct and error trials.');
end
pooledCorrect = sum(meanCorrect(includeMask,:) .* nCorrect(includeMask), 1) / totalCorrect;
pooledError = sum(meanError(includeMask,:) .* nError(includeMask), 1) / totalError;
kernelPC = pooledCorrect - pooledError;
end

% -------------------------------------------------------------------------
function assertUsableKernelSum(kernelSum, label)

if ~isfinite(kernelSum)
  error('makeBetaWeights:NonfiniteKernelSum', 'The %s kernel sum is nonfinite.', label);
end

if kernelSum <= 0
  error('makeBetaWeights:NonpositiveKernelSum', 'The %s kernel sum is %.6g; expected a positive value.', ...
    label, kernelSum);
end

% Guard against pathological normalization.
if abs(kernelSum) < 1e-9
  error('makeBetaWeights:NearZeroKernelSum', 'The %s kernel sum is too close to zero for stable normalization.', label);
end
end
