function fit = fitReadoutModelToScales(offsetsDeg, obsScale, obsVar, mt, modelName)
% fitReadoutModelToScales
%
% Fit mechanistic readout model to observed non-anchor scales.
%
% Inputs:
%   offsetsDeg : non-anchor offsets only, e.g. [45 180]
%   obsScale   : observed pooled scales at those offsets
%   obsVar     : bootstrap variances for those pooled scales
%   mt         : struct from makeMTKernel
%   modelName  : currently 'gaussian_offset' is the active default
%
% Output fields:
%   fit.modelName
%   fit.params
%   fit.offsetsFitDeg
%   fit.obsScale
%   fit.obsVar
%   fit.predScale
%   fit.resid
%   fit.weightedSSE
%   fit.exitFlag
%   fit.mtParams

offsetsDeg = offsetsDeg(:)';
obsScale   = obsScale(:)';
obsVar     = obsVar(:)';

assert(numel(offsetsDeg) == numel(obsScale), 'offsetsDeg and obsScale must match in length.');
assert(numel(offsetsDeg) == numel(obsVar),   'offsetsDeg and obsVar must match in length.');
assert(all(obsVar > 0), 'All obsVar entries must be positive.');

switch lower(modelName)
    case 'gaussian_offset'
        p0 = [30, 0];
        lb = [1, -2];
        ub = [180, 2];

    case 'dog'
        p0 = [20, 80, 0.3];
        lb = [1,  2,  0];
        ub = [180, 180, 5];

    otherwise
        error('fitReadoutModelToScales:UnknownModel', 'Unknown modelName: %s', modelName);
end

obj = @(p) weightedObjective(p, offsetsDeg, obsScale, obsVar, mt, modelName);

opts = optimset('Display', 'off');
[pHat, fval] = fminsearchbndLocal(obj, p0, lb, ub, opts);

predScale = predictScalesFromReadout(offsetsDeg, mt, modelName, pHat);
resid     = obsScale - predScale;

fit = struct();
fit.modelName    = modelName;
fit.params       = pHat;
fit.offsetsFitDeg = offsetsDeg;
fit.obsScale     = obsScale;
fit.obsVar       = obsVar;
fit.predScale    = predScale;
fit.resid        = resid;
fit.weightedSSE  = fval;
fit.exitFlag     = 1;
fit.mtParams     = struct( ...
    'sigmaDeg', mt.sigmaDeg, ...
    'nullRatioAbs', mt.nullRatioAbs, ...
    'phiDeg', mt.phiDeg);

end

function err = weightedObjective(params, offsetsDeg, obsScale, obsVar, mt, modelName)

pred = predictScalesFromReadout(offsetsDeg, mt, modelName, params);
resid = obsScale - pred;
err = sum((resid.^2) ./ obsVar);

if ~isfinite(err)
    err = 1e12;
end

% Gentle shape constraint for DOG: surround should be broader than center
if strcmpi(modelName, 'dog')
    sigmaC = params(1);
    sigmaS = params(2);
    if sigmaS <= sigmaC
        err = err + 1e6 + 1e3 * (sigmaC - sigmaS + 1);
    end
end

end

function [x,fval] = fminsearchbndLocal(fun, x0, LB, UB, options)

x0u = xtransformLocal(x0, LB, UB);
obj = @(xu) fun(invtransformLocal(xu, LB, UB));
[xu, fval] = fminsearch(obj, x0u, options);
x = invtransformLocal(xu, LB, UB);

end

function xu = xtransformLocal(x, LB, UB)
x  = x(:);
LB = LB(:);
UB = UB(:);
xu = zeros(size(x));

for i = 1:numel(x)
    if isfinite(LB(i)) && isfinite(UB(i))
        t = (x(i) - LB(i)) / (UB(i) - LB(i));
        t = min(max(t, 1e-8), 1 - 1e-8);
        xu(i) = log(t / (1 - t));
    elseif isfinite(LB(i))
        xu(i) = log(max(x(i) - LB(i), 1e-8));
    elseif isfinite(UB(i))
        xu(i) = log(max(UB(i) - x(i), 1e-8));
    else
        xu(i) = x(i);
    end
end
end

function x = invtransformLocal(xu, LB, UB)
xu = xu(:);
LB = LB(:);
UB = UB(:);
x = zeros(size(xu));

for i = 1:numel(xu)
    if isfinite(LB(i)) && isfinite(UB(i))
        t = 1 ./ (1 + exp(-xu(i)));
        x(i) = LB(i) + (UB(i) - LB(i)) * t;
    elseif isfinite(LB(i))
        x(i) = LB(i) + exp(xu(i));
    elseif isfinite(UB(i))
        x(i) = UB(i) - exp(xu(i));
    else
        x(i) = xu(i);
    end
end

x = x(:)';
end