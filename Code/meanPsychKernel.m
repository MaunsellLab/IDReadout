function [kernel, kernelVar] = meanPsychKernel(cohMat, trialOutcome, noiseAmpPC)
% computePsychophysicalKernel  Compute a psychophysical kernel over time.
%
%   [kernel, kernelSE] = computePsychophysicalKernel(cohMat, trialOutcome)
%
%   INPUTS:
%     cohMat       : m x nTrials matrix of coherence noise values
%                    (rows = time bins, columns = trials).
%     trialOutcome : nTrials x 1 vector of trial outcomes
%                    0 = correct, 1 = wrong.
%     noiseAmpPC   : amplitude of the binary coherence noise
%
%   OUTPUTS:
%     kernel   : m x 1 kernel, defined as
%                  mean(noise | wrong) - mean(noise | correct)
%     kernelSE : m x 1 standard error of the kernel at each time bin (based on known distribution variance)
%
%   Trials with trialOutcome not equal to 0 or 1 are ignored.

  % Basic input checks
  if size(cohMat, 2) ~= numel(trialOutcome)
    error('cohMat and trialOutcome must have the same number of trials (columns).');
  end
  
  trialOutcome = trialOutcome(:)'; % make row for logical indexing
  
  % Select valid trials: trialEnd must be 0 or 1
  valid = (trialOutcome == 0) | (trialOutcome == 1);
  if ~any(valid)
    error('computePsychophysicalKernel:NoValidTrials', ...
      'No trials with outcome 0 or 1 found.');
  end
  cohMat = cohMat(:, valid);
  trialOutcome = trialOutcome(valid);
  
  % Indices for correct and wrong trials
  idxCorrect = (trialOutcome == 0);
  idxWrong   = (trialOutcome == 1);
  nCorrect = sum(idxCorrect);
  nWrong   = sum(idxWrong);
  if nCorrect == 0 || nWrong == 0
    warning('computePsychophysicalKernel:Unbalanced', ...
      'One of the outcome classes (correct/wrong) is empty.');
  end
  
  % Means
  meanCorrect = mean(cohMat(:, idxCorrect), 2);
  meanWrong   = mean(cohMat(:, idxWrong),   2);
  kernel = meanCorrect - meanWrong;                 % correct minus wrong

  % Standard error of the kernel (per time bin)
  % Pref and probe can have different kernel var because it is affected by
  % the noise amplitude
  if nargout > 1
    sigma2 = noiseAmpPC^2;                          % Known stimulus variance for ±a, p = 0.5
    kernelVar = sigma2 * (1/nCorrect + 1/nWrong);   % Variance at each time bin
  end
end
