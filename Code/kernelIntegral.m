function   [kIntegrals, R, RVar] = kernelIntegral(kernels, kVars, msPerVFrame)

% integrate kernels
%
% Input:
%         kernels       2x2xm matrix of kernels 
%         kVars         2x2 matrix of kernel variance 
%         msPerVFrame   mm per video frame (kernel entry)
%
% Output:
%         kIntegrals    2x2 kernel integrals
%         R             1x2 vector of integral ratios
%         RVar          1x2 vector of integral ratio variances computed
%                         using variance-of-a-ratio approximation (delta method)
%
  [preStepMS, intStartMS, intDurMS] = integralWindowMS();
  kIntegrals = nan(2, 2);               % 2 x 2 kernel integral scalars (inc/dec, pref/probe)
  iIndices = round((preStepMS + intStartMS) / msPerVFrame) : round((preStepMS + intStartMS + intDurMS) / msPerVFrame);
  R = nan(1, 2);
  RVar = nan(1, 2);
  for p = 1:2
    for s = 1:2
      kIntegrals(s, p) = sum(kernels(s, p, iIndices)) / length(iIndices);
      R(s) = kIntegrals(s, 2) / kIntegrals(s, 1);
      RVar(s) = (kVars(s, 2) / (kIntegrals(s, 1)^2)) + (kIntegrals(s, 2)^2 * kVars(s, 1) / (kIntegrals(s, 1)^4));
    end
  end
