function fit = fitEffectiveWeightingToScales( ...
    offsetsDeg, obsScale, obsVar, mtModel, modelName, varargin)
% fitEffectiveWeightingToScales
%
% Fit a parameterized effective MT-to-choice weighting function to observed
% psychophysical scales at non-anchor offsets.
%
% Active model:
%   'dog'
%
% Baseline handling:
%   'BaselineMode','fixed'
%       fit params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain]
%       with baseline supplied separately via 'FixedBaseline'
%
%   'BaselineMode','fit'
%       fit params = [sigmaCenterDeg, sigmaSurroundDeg, surroundGain, baseline]
%
%   'BaselineMode','auto'
%       uses 'fixed' when the number of non-anchor offsets is smaller than
%       'MinOffsetsForFitBaseline', otherwise uses 'fit'

p = inputParser;
p.addParameter('BaselineMode', 'auto', @(x) ischar(x) || isstring(x));
p.addParameter('FixedBaseline', [], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
p.addParameter('MinOffsetsForFitBaseline', 3, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('Bounds', struct(), @isstruct);
p.parse(varargin{:});

opts = p.Results;

offsetsDeg = offsetsDeg(:)';
obsScale   = obsScale(:)';
obsVar     = obsVar(:)';

assert(numel(offsetsDeg) == numel(obsScale), ...
    'offsetsDeg and obsScale must match in length.');
assert(numel(offsetsDeg) == numel(obsVar), ...
    'offsetsDeg and obsVar must match in length.');
assert(all(obsVar > 0), ...
    'All obsVar entries must be positive.');

if ~strcmpi(modelName, 'dog')
    error('fitEffectiveWeightingToScales:UnsupportedModel', ...
        'Active pathway supports only modelName=''dog''.');
end

nNonAnchor = numel(offsetsDeg);
baselineMode = lower(char(string(opts.BaselineMode)));

switch baselineMode
    case 'auto'
        if nNonAnchor < opts.MinOffsetsForFitBaseline
            baselineModeResolved = 'fixed';
        else
            baselineModeResolved = 'fit';
        end

    case 'fixed'
        baselineModeResolved = 'fixed';

    case 'fit'
        baselineModeResolved = 'fit';

    otherwise
        error('fitEffectiveWeightingToScales:BadBaselineMode', ...
            'Unknown BaselineMode: %s', baselineMode);
end

switch baselineModeResolved
    case 'fixed'
        if isempty(opts.FixedBaseline)
            error('fitEffectiveWeightingToScales:MissingFixedBaseline', ...
                'BaselineMode=''fixed'' requires FixedBaseline.');
        end

        p0 = [20, 80, 0.3];
        lb = [1,  2,  0];
        ub = [180, 180, 5];

    case 'fit'
        p0 = [20, 80, 0.3, 0];
        lb = [1,  2,  0, -2];
        ub = [180, 180, 5,  2];
end

[lb, ub] = applyBoundsOverride(lb, ub, baselineModeResolved, opts.Bounds);

obj = @(p) weightedObjective( ...
    p, offsetsDeg, obsScale, obsVar, mtModel, ...
    modelName, baselineModeResolved, opts.FixedBaseline);

optimOpts = optimset('Display', 'off');
[pHat, fval] = fminsearchbndLocal(obj, p0, lb, ub, optimOpts);

switch baselineModeResolved
    case 'fixed'
        predScale = predictScalesFromEffectiveWeighting( ...
            offsetsDeg, mtModel, modelName, pHat, opts.FixedBaseline);
        fixedBaselineResolved = opts.FixedBaseline;

    case 'fit'
        predScale = predictScalesFromEffectiveWeighting( ...
            offsetsDeg, mtModel, modelName, pHat, []);
        fixedBaselineResolved = [];
end

resid = obsScale - predScale;

fit = struct();
fit.modelName      = modelName;
fit.params         = pHat;
fit.baselineMode   = baselineModeResolved;
fit.fixedBaseline  = fixedBaselineResolved;
fit.nFreeParams    = numel(pHat);
fit.offsetsFitDeg  = offsetsDeg;
fit.obsScale       = obsScale;
fit.obsVar         = obsVar;
fit.predScale      = predScale;
fit.resid          = resid;
fit.weightedSSE    = fval;
fit.exitFlag       = 1;
fit.mtForwardModelParams = struct( ...
    'sigmaDeg', mtModel.sigmaDeg, ...
    'nullRatioAbs', mtModel.nullRatioAbs, ...
    'phiDeg', mtModel.phiDeg);

end

function err = weightedObjective( ...
    params, offsetsDeg, obsScale, obsVar, mtModel, ...
    modelName, baselineModeResolved, fixedBaseline)

switch baselineModeResolved
    case 'fixed'
        pred = predictScalesFromEffectiveWeighting( ...
            offsetsDeg, mtModel, modelName, params, fixedBaseline);
        sigmaC = params(1);
        sigmaS = params(2);

    case 'fit'
        pred = predictScalesFromEffectiveWeighting( ...
            offsetsDeg, mtModel, modelName, params, []);
        sigmaC = params(1);
        sigmaS = params(2);
end

resid = obsScale - pred;
err = sum((resid.^2) ./ obsVar);

if ~isfinite(err)
    err = 1e12;
end

if sigmaS <= sigmaC
    err = err + 1e6 + 1e3 * (sigmaC - sigmaS + 1);
end

end

function [lb, ub] = applyBoundsOverride(lb, ub, baselineModeResolved, boundsStruct)

if ~isfield(boundsStruct, 'dog')
    return;
end

B = boundsStruct.dog;

switch baselineModeResolved
    case 'fixed'
        if isfield(B, 'lbFixed'), lb = B.lbFixed; end
        if isfield(B, 'ubFixed'), ub = B.ubFixed; end

    case 'fit'
        if isfield(B, 'lbFit'), lb = B.lbFit; end
        if isfield(B, 'ubFit'), ub = B.ubFit; end
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