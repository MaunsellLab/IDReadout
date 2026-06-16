function deltaM = mtReadoutTemplate(deltaDeg, mtModel, varargin)
% mtReadoutTemplate  Effective MT perturbation for one probe offset.
%
% For 0 < delta < 180 deg, the stored probe predictor represents the sum of
% two equally yoked streams. Per unit of that summed predictor, the effective
% MT template is the average of the +delta and -delta single-stream effects.
%
% Name-value option:
%   'TemplateMode'  'signed' (default) or 'rectified'

p = inputParser;
addParameter(p, 'TemplateMode', 'signed', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
templateMode = lower(char(string(p.Results.TemplateMode)));

phiDeg = mtModel.phiDeg(:)';
deltaDeg = abs(double(deltaDeg));

switch templateMode
    case 'signed'
        if isSingleStreamOffset(deltaDeg)
            deltaM = computeMTDeltaM(phiDeg, deltaDeg, mtModel);
        else
            plusTemplate  = computeMTDeltaM(phiDeg, +deltaDeg, mtModel);
            minusTemplate = computeMTDeltaM(phiDeg, -deltaDeg, mtModel);
            deltaM = 0.5 .* (plusTemplate + minusTemplate);
        end

    case 'rectified'
        if isSingleStreamOffset(deltaDeg)
            deltaM = max(computeMTDeltaM(phiDeg, deltaDeg, mtModel), 0);
        else
            plusTemplate  = max(computeMTDeltaM(phiDeg, +deltaDeg, mtModel), 0);
            minusTemplate = max(computeMTDeltaM(phiDeg, -deltaDeg, mtModel), 0);
            deltaM = 0.5 .* (plusTemplate + minusTemplate);
        end

    otherwise
        error('mtReadoutTemplate:UnknownMode', ...
            'Unknown MT template mode: %s', templateMode);
end
end

function tf = isSingleStreamOffset(deltaDeg)
tf = abs(deltaDeg) < 1e-9 || abs(deltaDeg - 180) < 1e-9;
end

function deltaM = computeMTDeltaM(phiDeg, deltaDeg, mtModel)
if numel(phiDeg) ~= numel(mtModel.phiDeg) || any(phiDeg ~= mtModel.phiDeg)
    error('mtReadoutTemplate:GridMismatch', ...
        'The supplied MT model must use one shared phi grid.');
end

shifted = mod(phiDeg - deltaDeg + 180, 360) - 180;
G = exp(-(shifted.^2) ./ (2 * mtModel.sigmaMTDeg.^2));
deltaM = G - mean(G);
end
