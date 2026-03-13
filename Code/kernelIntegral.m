function [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame)

  [preStepMS, intStartMS, intDurMS] = integralWindowMS();
  
  % 5 x 2 x 2 kernel integrals (diff/c/nC/RF/Opp, inc/dec, pref/probe)
  kIntegrals = nan(5, 2, 2);
  iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
    round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
  N = numel(iIndices);
  
  % 5 x 2 statistics (diff/c/nC/RF/Opp, inc/dec)
  R    = nan(5, 2);
  RVar = nan(5, 2);
  
  for sideType = 1:5
    for s = 1:2
      for p = 1:2
        kIntegrals(sideType, s, p) = sum(kernels(sideType, s, p, iIndices)) / N;
      end
  
      % Approx variance of the mean over N bins
      varPrefInt  = kVars(s, 1) / N;
      varProbeInt = kVars(s, 2) / N;
  
      % Ratio and delta-method variance
      denom = kIntegrals(sideType, s, 1);
      R(sideType, s) = kIntegrals(sideType, s, 2) / denom;
  
      RVar(sideType, s) = (varProbeInt / (denom^2)) + ...
        ((kIntegrals(sideType, s, 2)^2) * varPrefInt / (denom^4));
    end
  end
end