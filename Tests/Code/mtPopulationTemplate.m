function deltaM = mtPopulationTemplate(directionDeg,mtModel)
% mtPopulationTemplate  Signed MT perturbation for individual directions.
%
% Rows of deltaM correspond to directionDeg. Each row is the response to
% one unit of physical coherence at that direction:
%   G(phi-direction) - mean_phi G(phi-direction)
%
% Unlike mtReadoutTemplate, this function never averages yoked probes.

required = {'phiDeg','sigmaMTDeg'};
if ~isstruct(mtModel) || ~all(isfield(mtModel,required))
  error('mtPopulationTemplate:BadMTModel', ...
    'mtModel must contain phiDeg and sigmaMTDeg.');
end
phiDeg = double(mtModel.phiDeg(:)');
directionDeg = double(directionDeg(:));
if any(~isfinite(directionDeg)) || any(~isfinite(phiDeg)) || ...
    ~isscalar(mtModel.sigmaMTDeg) || ~isfinite(mtModel.sigmaMTDeg) || ...
    mtModel.sigmaMTDeg<=0
  error('mtPopulationTemplate:InvalidInput','Directions and MT model must be finite.');
end

distanceDeg = circularDifference(phiDeg,directionDeg);
G = exp(-(distanceDeg.^2)./(2*double(mtModel.sigmaMTDeg).^2));
deltaM = G-mean(G,2);
end

function d = circularDifference(phiDeg,directionDeg)
d = mod(phiDeg-directionDeg+180,360)-180;
end
