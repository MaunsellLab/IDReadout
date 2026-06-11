function wEff = makeEffectiveWeighting(phiDeg, modelName, params, fixedBaseline)
% makeEffectiveWeighting
%
% Construct the effective MT-to-choice weighting function wEff(phi).
%
% Supported models:
%   'dog'
%       params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain]
%           uses fixedBaseline as the asymptotic baseline
%
%       params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain, baseline]
%           baseline is taken from params(4), and fixedBaseline is ignored
%
% Notes:
%   - Overall gain is irrelevant because predicted scales are normalized by
%     the value at 0 deg in the forward model.
%   - The DOG family is interpreted as a parameterized family for the
%     effective weighting function itself.

if nargin < 4
    fixedBaseline = [];
end

phiDeg = wrapTo180Local(phiDeg);

switch lower(modelName)
    case 'dog'
        validateattributes(params, {'numeric'}, {'vector','real'});
        nParams = numel(params);

        if nParams == 3
            sigmaC = params(1);
            sigmaS = params(2);
            alpha  = params(3);

            if isempty(fixedBaseline)
                error('makeEffectiveWeighting:MissingFixedBaseline', ...
                    ['DOG with 3 parameters requires fixedBaseline. ' ...
                     'Supply fixedBaseline or pass a 4-parameter DOG.']);
            end
            baseline = fixedBaseline;

        elseif nParams == 4
            sigmaC   = params(1);
            sigmaS   = params(2);
            alpha    = params(3);
            baseline = params(4);

        else
            error('makeEffectiveWeighting:BadDOGParamCount', ...
                'DOG requires 3 or 4 parameters.');
        end

        hRaw = exp(-(phiDeg.^2) ./ (2 * sigmaC^2)) ...
             - alpha .* exp(-(phiDeg.^2) ./ (2 * sigmaS^2));

        h0 = 1 - alpha;
        if ~isfinite(h0) || abs(h0) < 1e-9
            error('makeEffectiveWeighting:BadDOGNormalization', ...
                'DOG normalization at 0 deg is invalid.');
        end

        h = hRaw ./ h0;
        wEff = baseline + (1 - baseline) .* h;

    otherwise
        error('makeEffectiveWeighting:UnknownModel', ...
            'Unknown modelName: %s', modelName);
end

end

function x = wrapTo180Local(x)
x = mod(x + 180, 360) - 180;
end