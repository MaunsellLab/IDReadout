function [H, meanT] = centeredMTSensitivity(deltaDeg, sigmaDeg, centerTuning)
% CENTEREDMTSENSITIVITY  Wrapped Gaussian MT sensitivity, optionally centered.
%
% deltaDeg is wrapped to [-180,180).  The circular mean is the exact integral
% of the wrapped-distance Gaussian over that interval divided by 360 deg.

deltaDeg = mod(deltaDeg + 180, 360) - 180;
T = exp(-0.5*(deltaDeg./sigmaDeg).^2);
meanT = sigmaDeg*sqrt(2*pi)/360 * erf(180/(sqrt(2)*sigmaDeg));
if centerTuning
    H = T - meanT;
else
    H = T;
end
end

