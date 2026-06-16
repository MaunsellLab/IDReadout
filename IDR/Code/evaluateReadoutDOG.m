function [aPhi, paramStruct] = evaluateReadoutDOG(phiDeg, params)
% evaluateReadoutDOG  Evaluate a DOG readout over MT preferred direction.
%
% params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain]
% or, for the rectified-template model,
% params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain, baselineOffset]

phiDeg = phiDeg(:)';
params = params(:)';

if numel(params) == 3
    baselineOffset = 0;
elseif numel(params) == 4
    baselineOffset = params(4);
else
    error('evaluateReadoutDOG:BadParamCount', ...
        'Expected three or four DOG parameters, got %d.', numel(params));
end

sigmaC = params(1);
sigmaS = params(2);
As = params(3);

aPhi = baselineOffset ...
    + exp(-(phiDeg.^2) ./ (2 * sigmaC.^2)) ...
    - As .* exp(-(phiDeg.^2) ./ (2 * sigmaS.^2));

paramStruct = struct( ...
    'sigmaCenterDeg', sigmaC, ...
    'sigmaSurroundDeg', sigmaS, ...
    'surroundGain', As, ...
    'baselineOffset', baselineOffset);
end
