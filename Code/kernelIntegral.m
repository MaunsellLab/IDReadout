function [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame)

  [preStepMS, intStartMS, intDurMS] = integralWindowMS();

  nSideTypes = size(kernels, 1);

  % nSideTypes x 2 x 2 kernel integrals (sideType, inc/dec, pref/probe)
  kIntegrals = nan(nSideTypes, 2, 2);

  iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
    round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
  N = numel(iIndices);

  % nSideTypes x 2 statistics (sideType, inc/dec)
  R    = nan(nSideTypes, 2);
  RVar = nan(nSideTypes, 2);

  for sideType = 1:nSideTypes
    for s = 1:2
      for p = 1:2
        kIntegrals(sideType, s, p) = ...
          sum(kernels(sideType, s, p, iIndices)) / N;
      end

      % Approx variance of the mean over N bins
      varPrefInt  = kVars(sideType, s, 1) / N;
      varProbeInt = kVars(sideType, s, 2) / N;

      % Ratio and delta-method variance
      denom = kIntegrals(sideType, s, 1);

      if ~isfinite(denom) || denom == 0
        continue
      end

      R(sideType, s) = kIntegrals(sideType, s, 2) / denom;

      RVar(sideType, s) = (varProbeInt / (denom^2)) + ...
        ((kIntegrals(sideType, s, 2)^2) * varPrefInt / (denom^4));
    end
  end
end