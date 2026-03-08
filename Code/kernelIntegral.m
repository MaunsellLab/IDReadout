function [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame)

  [preStepMS, intStartMS, intDurMS] = integralWindowMS();
  kIntegrals = nan(2, 2);

  iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
             round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
  N = numel(iIndices);

  R    = nan(1, 2);
  RVar = nan(1, 2);

  for s = 1:2
    for p = 1:2
      kIntegrals(s, p) = sum(kernels(s, p, iIndices)) / N;
    end

    % Approx variance of the mean over N bins
    varPrefInt  = kVars(s, 1) / N;
    varProbeInt = kVars(s, 2) / N;

    % Ratio and delta-method variance
    denom = kIntegrals(s, 1);
    R(s) = kIntegrals(s, 2) / denom;

    RVar(s) = (varProbeInt / (denom^2)) + ...
              ((kIntegrals(s, 2)^2) * varPrefInt / (denom^4));
  end
end