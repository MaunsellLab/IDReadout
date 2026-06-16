function mtModel = makeMTReadoutForwardModel(varargin)
% makeMTReadoutForwardModel  Construct the fixed MT population-response model.
%
% The same phi grid is used for MT tuning, readout evaluation, and discrete
% mean subtraction. The grid should not contain both -180 and +180 degrees.

p = inputParser;
addParameter(p, 'sigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'phiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
parse(p, varargin{:});

phiDeg = p.Results.phiDeg(:)';
sigmaMTDeg = p.Results.sigmaMTDeg;

mtModel = struct();
mtModel.sigmaMTDeg = sigmaMTDeg;
mtModel.phiDeg = phiDeg;
mtModel.gridStepDeg = median(diff(phiDeg));

G0 = exp(-(phiDeg.^2) ./ (2 * sigmaMTDeg.^2));
mtModel.G0 = G0;
mtModel.Gbar = mean(G0);
mtModel.note = ['Delta m(phi;delta) is defined as G(phi-delta) minus the ' ...
    'discrete mean of G on the same phi-grid.'];
end
