function [kernel, kernelVar, stats] = meanPsychKernel(cohMat, trialOutcome, changeSide, noiseAmpPC)
% meanPsychKernel  Compute a psychophysical kernel over time.
%
%   [kernel, kernelVar] = meanPsychKernel(cohMat, trialOutcome, noiseAmpPC)
%   [kernel, kernelVar, stats] = meanPsychKernel(...)
%
%   INPUTS:
%     cohMat       : m x nTrials matrix of coherence noise values
%                    (rows = time bins, columns = trials).
%     trialOutcome : nTrials x 1 (or 1 x nTrials) vector of trial outcomes
%                    0 = correct (hit), 1 = wrong (miss).
%     noiseAmpPC   : amplitude of the binary coherence noise (in percent coherence units)
%
%   OUTPUTS:
%     kernel    : m x 1 kernel, defined as mean(correct) - mean(wrong)
%     kernelVar : scalar variance of the kernel at each time bin (known stimulus variance)
%     stats     : (optional) sufficient statistics to allow grand-pooled kernels:
%                .nCorrect, .nWrong, .sumCorrect (m x 1), .sumWrong (m x 1)
%
%   Trials with trialOutcome not equal to 0 or 1 are ignored.

  % ---- Basic input checks ----
  if size(cohMat, 2) ~= numel(trialOutcome)
    error('cohMat and trialOutcome must have the same number of trials (columns).');
  end

  trialOutcome = trialOutcome(:)'; % row for logical indexing

  % ---- Select valid trials ----
  valid = (trialOutcome == 0) | (trialOutcome == 1);
  if ~any(valid)
    error('meanPsychKernel:NoValidTrials', ...
      'No trials with outcome 0 or 1 found.');
  end
  cohMat = cohMat(:, valid);
  trialOutcome = trialOutcome(valid);

  % ---- Indices for correct and wrong trials ----
  idxCorrect = (trialOutcome == 0);
  idxWrong   = (trialOutcome == 1);
  nCorrect = sum(idxCorrect);
  nWrong   = sum(idxWrong);

  if nCorrect == 0 || nWrong == 0
    warning('meanPsychKernel:Unbalanced', ...
      'One of the outcome classes (correct/wrong) is empty.');
  end

  % ---- Sufficient statistics and means ----
  sumCorrect = sum(cohMat(:, idxCorrect), 2);
  sumWrong   = sum(cohMat(:, idxWrong),   2);

  if nCorrect > 0
    meanCorrect = sumCorrect ./ nCorrect;
  else
    meanCorrect = nan(size(sumCorrect));
  end
  if nWrong > 0
    meanWrong = sumWrong ./ nWrong;
  else
    meanWrong = nan(size(sumWrong));
  end

  kernel = meanCorrect - meanWrong;    % correct minus wrong (matches existing convention)

  % ---- Variance (scalar; same for every time bin under iid ±a noise) ----
  sigma2 = noiseAmpPC^2;                          % Known stimulus variance for ±a, p = 0.5
  kernelVar = sigma2 * (1/max(nCorrect,1) + 1/max(nWrong,1));

  % ---- Optional output for grand pooling ----
  if nargout > 2
    stats = struct();
    idxRF   = (changeSide == 0);
    stats.nRFCorrect  = sum(trialOutcome(idxRF) == 0);
    stats.nRFWrong    = sum(trialOutcome(idxRF) ~= 0);
    stats.nCorrect    = nCorrect;
    stats.nWrong      = nWrong;
    stats.sumCorrect  = sumCorrect;
    stats.sumWrong    = sumWrong;
    stats.sigma2      = sigma2;

  end
end
