function deltaM = mtPopulationTemplate(directionDeg,mtModel,varargin)
% mtPopulationTemplate  Signed MT perturbation for individual directions.
%
% Rows of deltaM correspond to directionDeg. Each row is the response to
% one unit of physical coherence at that direction:
%   G(phi-direction) - mean_phi G(phi-direction)
%
% Unlike mtReadoutTemplate, this function never averages yoked probes.
%
% Name-value option:
%   MTTemplateMode  "signed" (default) or "rectified". Rectified mode
%   applies max(deltaM,0) to each unit positive-coherence direction
%   template before its signed physical coherence coefficient is applied by
%   downstream code. It is not an absolute firing-rate floor.

p=inputParser;
addParameter(p,'MTTemplateMode',"signed",@validTemplateMode);
parse(p,varargin{:});mode=lower(string(p.Results.MTTemplateMode));

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
if mode=="rectified",deltaM=max(deltaM,0);end
end

function d = circularDifference(phiDeg,directionDeg)
d = mod(phiDeg-directionDeg+180,360)-180;
end

function tf=validTemplateMode(x)
tf=(ischar(x)||isstring(x))&&isscalar(string(x))&& ...
  any(strcmpi(string(x),["signed","rectified"]));
end
