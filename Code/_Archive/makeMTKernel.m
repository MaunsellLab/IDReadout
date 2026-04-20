function mt = makeMTKernel(varargin)
% makeMTKernel
%
% Build fixed MT directional response-change kernel:
%   g(theta) = A * exp(-theta^2 / (2*sigma^2)) + B
%
% Constraints:
%   g(0)   = 1
%   g(180) = -nullRatioAbs
%
% Default assumptions:
%   sigmaDeg     = 37.5
%   nullRatioAbs = 1/3
%   phiDeg       = -180:1:179
%
% Example:
%   mt = makeMTKernel();
%
% Returns struct:
%   mt.sigmaDeg
%   mt.nullRatioAbs
%   mt.phiDeg
%   mt.A
%   mt.B
%   mt.gFun              % function handle g(thetaDeg)
%   mt.shiftedGFun       % function handle m(phiDeg, deltaDeg)

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

mt = struct();
mt.sigmaDeg     = sigmaDeg;
mt.nullRatioAbs = nullRatioAbs;
mt.phiDeg       = phiDeg;
mt.A            = A;
mt.B            = B;
mt.gFun         = @(thetaDeg) mtKernelEval(thetaDeg, sigmaDeg, A, B);
mt.shiftedGFun  = @(phiDegLocal, deltaDeg) mtKernelEval(wrapTo180Local(phiDegLocal - deltaDeg), sigmaDeg, A, B);

end

function g = mtKernelEval(thetaDeg, sigmaDeg, A, B)
thetaDeg = wrapTo180Local(thetaDeg);
g = A .* exp(-(thetaDeg.^2) ./ (2 * sigmaDeg^2)) + B;
end

function x = wrapTo180Local(x)
x = mod(x + 180, 360) - 180;
end