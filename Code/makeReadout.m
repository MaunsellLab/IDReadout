function w = makeReadout(phiDeg, modelName, params)
% makeReadout
%
% Construct latent readout function w(phi).
%
% Supported models:
%   'gaussian_offset' : params = [sigmaW, offset]
%   'dog'             : params = [sigmaCenter, sigmaSurround, alpha]
%                       w = exp(-phi^2/(2*s1^2)) - alpha*exp(-phi^2/(2*s2^2))
%
% Notes:
%   - Overall gain is irrelevant because predicted scales are normalized by R(0).
%   - DOG is scaffolded for future use but need not be fit yet.

phiDeg = wrapTo180Local(phiDeg);

switch lower(modelName)
    case 'gaussian_offset'
        validateattributes(params, {'numeric'}, {'vector','numel',2,'real'});
        sigmaW = params(1);
        c      = params(2);

        w = exp(-(phiDeg.^2) ./ (2 * sigmaW^2)) + c;

    case 'dog'
        validateattributes(params, {'numeric'}, {'vector','numel',3,'real'});
        sigmaC = params(1);
        sigmaS = params(2);
        alpha  = params(3);

        w = exp(-(phiDeg.^2) ./ (2 * sigmaC^2)) ...
          - alpha .* exp(-(phiDeg.^2) ./ (2 * sigmaS^2));

    otherwise
        error('makeReadout:UnknownModel', 'Unknown modelName: %s', modelName);
end

end

function x = wrapTo180Local(x)
x = mod(x + 180, 360) - 180;
end