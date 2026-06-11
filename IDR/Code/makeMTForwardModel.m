function mtModel = makeMTForwardModel(varargin)
% makeMTForwardModel
%
% Build a fixed MT forward model for directional response change across MT
% preferred-direction channels.
%
% Forward model:
%   rMT(theta) = A * exp(-theta^2 / (2*sigma^2)) + B
%
% Constraints:
%   rMT(0)   = 1
%   rMT(180) = -nullRatioAbs
%
% Default assumptions:
%   sigmaDeg     = 37.5
%   nullRatioAbs = 1/3
%   phiDeg       = -180:1:179
%
% Example:
%   mtModel = makeMTForwardModel();
%
% Returns struct:
%   mtModel.sigmaDeg
%   mtModel.nullRatioAbs
%   mtModel.phiDeg
%   mtModel.A
%   mtModel.B
%   mtModel.responseFun         % function handle rMT(thetaDeg)
%   mtModel.shiftedResponseFun  % function handle rMT(phiDeg, deltaDeg)
%
p = inputParser;
p.addParameter('sigmaDeg', 37.5, @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
p.addParameter('nullRatioAbs', 1/3, @(x) validateattributes(x, {'numeric'}, {'scalar','real','positive'}));
p.addParameter('phiDeg', -180:1:179, @(x) validateattributes(x, {'numeric'}, {'vector','real'}));
p.parse(varargin{:});

sigmaDeg     = p.Results.sigmaDeg;
nullRatioAbs = p.Results.nullRatioAbs;
phiDeg       = p.Results.phiDeg(:)';

e180 = exp(-(180^2) / (2 * sigmaDeg^2));

% Solve:
%   A + B      = 1
%   A*e180 + B = -nullRatioAbs
A = (1 + nullRatioAbs) / (1 - e180);
B = 1 - A;

mtModel = struct();
mtModel.sigmaDeg     = sigmaDeg;
mtModel.nullRatioAbs = nullRatioAbs;
mtModel.phiDeg       = phiDeg;
mtModel.A            = A;
mtModel.B            = B;

mtModel.responseFun        = @(thetaDeg) mtForwardEval(thetaDeg, sigmaDeg, A, B);
mtModel.shiftedResponseFun = @(phiDegLocal, deltaDeg) ...
    mtForwardEval(wrapTo180Local(phiDegLocal - deltaDeg), sigmaDeg, A, B);
end

function r = mtForwardEval(thetaDeg, sigmaDeg, A, B)
thetaDeg = wrapTo180Local(thetaDeg);
r = A .* exp(-(thetaDeg.^2) ./ (2 * sigmaDeg^2)) + B;
end

function x = wrapTo180Local(x)
x = mod(x + 180, 360) - 180;
end