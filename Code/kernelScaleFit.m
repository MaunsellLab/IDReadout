function [scale, scaleSEM, fitR2, sse] = kernelScaleFit(kernels, msPerVFrame)
% Fit probe kernel as a scaled version of pref kernel over the integration window.
%
% Model:
%       probe(t) ~= scale * pref(t)
%
% Uses a through-origin least-squares fit over the same time bins used by
% kernelIntegral():
%
%       scale = sum(pref .* probe) / sum(pref.^2)
%
% Input:
%       kernels        nSideTypes x 2 x 2 x T array
%                      (sideType, inc/dec, pref/probe, time)
%       msPerVFrame    ms per video frame (kernel entry)
%
% Output:
%       scale          nSideTypes x 2 matrix of fitted scale factors
%       scaleSEM       nSideTypes x 2 matrix of approximate SEMs for scale
%       fitR2          nSideTypes x 2 matrix of through-origin R^2 values
%       sse            nSideTypes x 2 matrix of residual sums of squares

  [preStepMS, intStartMS, intDurMS] = integralWindowMS();

  iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
    round((preStepMS + intStartMS + intDurMS) / msPerVFrame);

  nBins = numel(iIndices);
  nSideTypes = size(kernels, 1);

  scale    = nan(nSideTypes, 2);
  scaleSEM = nan(nSideTypes, 2);
  fitR2    = nan(nSideTypes, 2);
  sse      = nan(nSideTypes, 2);

  for sideType = 1:nSideTypes
    for s = 1:2

      pref  = squeeze(kernels(sideType, s, 1, iIndices));
      probe = squeeze(kernels(sideType, s, 2, iIndices));

      % Force column vectors
      pref  = pref(:);
      probe = probe(:);

      xx = sum(pref .^ 2);

      if xx == 0 || any(~isfinite(pref)) || any(~isfinite(probe))
        continue
      end

      % Through-origin least-squares scale factor
      scale(sideType, s) = sum(pref .* probe) / xx;

      % Residuals and SSE
      resid = probe - scale(sideType, s) * pref;
      sse(sideType, s) = sum(resid .^ 2);

      % Approximate residual variance for through-origin regression
      % df = nBins - 1 because one parameter (scale) is estimated
      if nBins > 1
        sigma2 = sse(sideType, s) / (nBins - 1);
        scaleSEM(sideType, s) = sqrt(sigma2 / xx);
      end

      % Through-origin R^2:
      % fraction of probe sum-of-squares explained by scaled pref
      denom = sum(probe .^ 2);
      if denom > 0
        fitR2(sideType, s) = 1 - sse(sideType, s) / denom;
      end
    end
  end
end