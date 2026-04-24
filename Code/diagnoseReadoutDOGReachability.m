% ========================================================================
function diag = diagnoseReadoutDOGReachability(varargin)
% Brute-force diagnostic for the DOG readout model.
%
% Maps the reachable (S45, S180) pairs over a parameter grid to determine
% whether the empirical target lies within the model family's range.
%
% Example:
%   diag = diagnoseReadoutDOGReachability( ...
%       'TargetS45', 0.206, ...
%       'TargetS180', 0.584, ...
%       'MakePlot', true);

p = inputParser;
addParameter(p, 'SigmaCenterDeg', 5:5:120, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SigmaSurroundDeg', 10:5:180, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SurroundGain', 0:0.05:2.0, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'SigmaMTDeg', 37.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'PhiDeg', -180:1:179, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'TargetS45', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'TargetS180', NaN, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MakePlot', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opt = p.Results;

mtModel = makeMTReadoutForwardModel( ...
    'sigmaMTDeg', opt.SigmaMTDeg, ...
    'phiDeg', opt.PhiDeg);

sC = opt.SigmaCenterDeg(:)';
sS = opt.SigmaSurroundDeg(:)';
aS = opt.SurroundGain(:)';

nTotal = numel(sC) * numel(sS) * numel(aS);

S45 = nan(nTotal, 1);
S180 = nan(nTotal, 1);
paramsMat = nan(nTotal, 3);
validMask = false(nTotal, 1);

idx = 0;
for i = 1:numel(sC)
    for j = 1:numel(sS)
        if sS(j) <= sC(i)
            continue;
        end
        for k = 1:numel(aS)
            idx = idx + 1;
            params = [sC(i), sS(j), aS(k)];
            pred = predictNormalizedScaleFromReadout([45 180], mtModel, params);

            paramsMat(idx, :) = params;
            S45(idx) = pred(1);
            S180(idx) = pred(2);
            validMask(idx) = all(isfinite(pred));
        end
    end
end

paramsMat = paramsMat(1:idx, :);
S45 = S45(1:idx);
S180 = S180(1:idx);
validMask = validMask(1:idx);

paramsMat = paramsMat(validMask, :);
S45 = S45(validMask);
S180 = S180(validMask);

diag = struct();
diag.paramNames = {'sigmaCenterDeg', 'sigmaSurroundDeg', 'surroundGain'};
diag.params = paramsMat;
diag.S45 = S45;
diag.S180 = S180;
diag.target = [opt.TargetS45, opt.TargetS180];

if isfinite(opt.TargetS45) && isfinite(opt.TargetS180) && ~isempty(S45)
    d2 = (S45 - opt.TargetS45).^2 + (S180 - opt.TargetS180).^2;
    [bestD2, bestIdx] = min(d2);
    diag.bestIdx = bestIdx;
    diag.bestParams = paramsMat(bestIdx, :);
    diag.bestPrediction = [S45(bestIdx), S180(bestIdx)];
    diag.bestDistance = sqrt(bestD2);
else
    diag.bestIdx = [];
    diag.bestParams = [];
    diag.bestPrediction = [];
    diag.bestDistance = NaN;
end

if opt.MakePlot && ~isempty(S45)
    fig = figure; clf; hold on;
    scatter(S45, S180, 10, 'filled');
    xlabel('Predicted S(45)');
    ylabel('Predicted S(180)');
    title('Reachable (S45, S180) pairs for DOG readout model');
    box off;

    if isfinite(opt.TargetS45) && isfinite(opt.TargetS180)
        plot(opt.TargetS45, opt.TargetS180, 'kp', 'MarkerSize', 14, ...
            'MarkerFaceColor', 'k');
    end
end
end