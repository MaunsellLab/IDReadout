function F = predictKernelTrials(beta1, pCorrect, stdSumDM, prefAmp, probeAmp, proj, msPerVFrame, varargin)
% predictKernelTrials  Predict trials needed for preferred/probe kernels
% to exceed zero by a one-sided z-threshold, using a scale-invariant forecast
% based on session m1 and an assumed active integration window.
%
% Required inputs:
%   beta1        slope from m1: isCorrect ~ sumDM
%   pCorrect     session mean P(correct)
%   stdSumDM     SD of the session's sumDM (the actual model predictor)
%   prefAmp      binary preferred-stream amplitude (e.g. 5 for ±5% coherence)
%   probeAmp     binary probe-stream amplitude
%   proj         probe projection factor used in sumDM construction
%   msPerVFrame  ms per video frame
%
% Optional name/value:
%   'refWindowMS'  reference window for m1 / sumDM, default 250
%   'windowMS'     active window(s) for forecast, default [50 125 250]
%   'zTarget'      one-sided threshold in SEM units, default 2
%   'maxN'         cap for predicted N, default 1e9
%
% Output struct F:
%   F.windowMS
%   F.nFrames
%   F.prefN
%   F.probeN
%   F.gamma
%   F.rPref
%   F.rProbe
%   F.valid
%
% Notes:
%   gamma = beta1 * stdSumDM is the scale-invariant sensitivity.
%   The forecast assumes the kernel area implied by m1 is spread uniformly
%   over the active window. Shorter active windows therefore predict larger
%   kernel heights and smaller required N.

  p = inputParser;
  p.addRequired('beta1', @(x) isnumeric(x) && isscalar(x));
  p.addRequired('pCorrect', @(x) isnumeric(x) && isscalar(x));
  p.addRequired('stdSumDM', @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.addRequired('prefAmp', @(x) isnumeric(x) && isscalar(x) && x >= 0);
  p.addRequired('probeAmp', @(x) isnumeric(x) && isscalar(x) && x >= 0);
  p.addRequired('proj', @(x) isnumeric(x) && isscalar(x));
  p.addRequired('msPerVFrame', @(x) isnumeric(x) && isscalar(x) && x > 0);

  p.addParameter('refWindowMS', 250, @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.addParameter('windowMS', [50 125 250], @(x) isnumeric(x) && isvector(x) && all(x > 0));
  p.addParameter('zTarget', 2, @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.addParameter('maxN', 1e9, @(x) isnumeric(x) && isscalar(x) && x > 0);
  p.parse(beta1, pCorrect, stdSumDM, prefAmp, probeAmp, proj, msPerVFrame, varargin{:});

  refWindowMS = p.Results.refWindowMS;
  windowMS = p.Results.windowMS(:)';
  zTarget = p.Results.zTarget;
  maxN = p.Results.maxN;

  % Scale-invariant standardized sensitivity
  gamma = beta1 * stdSumDM;

  % Relative contribution of preferred / probe to the combined decision axis
  sigmaComb = sqrt(prefAmp^2 + (proj * probeAmp)^2);
  if sigmaComb <= 0 || ~isfinite(sigmaComb)
    rPref = NaN;
    rProbe = NaN;
  else
    rPref  = prefAmp / sigmaComb;
    rProbe = abs(proj) * probeAmp / sigmaComb;
  end

  nRefFrames = max(1, round(refWindowMS / msPerVFrame));
  nFrames = max(1, round(windowMS / msPerVFrame));

  F = struct();
  F.windowMS = windowMS;
  F.nFrames = nFrames;
  F.refWindowMS = refWindowMS;
  F.nRefFrames = nRefFrames;
  F.gamma = gamma;
  F.rPref = rPref;
  F.rProbe = rProbe;
  F.beta1 = beta1;
  F.pCorrect = pCorrect;
  F.stdSumDM = stdSumDM;
  F.zTarget = zTarget;
  F.valid = true;

  if ~isfinite(gamma) || gamma <= 0 || ~isfinite(pCorrect) || pCorrect <= 0 || pCorrect >= 1
    F.prefN = nan(size(windowMS));
    F.probeN = nan(size(windowMS));
    F.valid = false;
    return;
  end

  % This is the conservative "kernel height" forecast:
  % N_target = [ z * T_active / (gamma * r_j * sqrt(T_ref) * [p(1-p)]^(3/2)) ]^2
  commonDen = gamma * sqrt(nRefFrames) * (pCorrect * (1 - pCorrect))^(3/2);

  F.prefN = localPredictHeightN(zTarget, nFrames, commonDen, rPref, maxN);
  F.probeN = localPredictHeightN(zTarget, nFrames, commonDen, rProbe, maxN);

end


function N = localPredictHeightN(zTarget, nFrames, commonDen, rj, maxN)

  if ~isfinite(rj) || rj <= 0
    N = nan(size(nFrames));
    return;
  end

  denom = commonDen * rj;
  N = (zTarget .* nFrames ./ denom).^2;

  bad = ~isfinite(N) | N < 0;
  N(bad) = NaN;
  N(N > maxN) = maxN;
end