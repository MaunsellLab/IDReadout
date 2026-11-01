function G = circularTuning(mechanismDeg, streamDeg, sigmaDeg)
%CIRCULARTUNING Circular Gaussian direction tuning.
delta = mod(mechanismDeg(:) - streamDeg(:)' + 180, 360) - 180;
G = exp(-(delta.^2) ./ (2*sigmaDeg.^2));
end

