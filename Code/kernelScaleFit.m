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
%       kernels        5 x 2 x 2 x T array
%                      (diff/c/nC/RF/Opp, inc/dec, pref/probe, time)
%       msPerVFrame    ms per video frame (kernel entry)
%
% Output:
%       scale          5 x 2 matrix of fitted scale factors
%       scaleSEM       5 x 2 matrix of approximate SEMs for scale
%       fitR2          5 x 2 matrix of through-origin R^2 values
%       sse            5 x 2 matrix of residual sums of squares
%
% Notes:
%   1. This treats pref as the predictor and probe as the response.
%   2. SEM is based on standard OLS residual variance, assuming independent,
%      homoscedastic residuals across time bins and negligible error in pref.
%   3. Because both kernels are themselves estimated quantities, this SEM is
%      only approximate. It is still a useful descriptive measure.

[preStepMS, intStartMS, intDurMS] = integralWindowMS();

iIndices = round((preStepMS + intStartMS) / msPerVFrame) : ...
  round((preStepMS + intStartMS + intDurMS) / msPerVFrame);

nBins = numel(iIndices);

scale    = nan(5, 2);
scaleSEM = nan(5, 2);
fitR2    = nan(5, 2);
sse      = nan(5, 2);

for sideType = 1:5
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