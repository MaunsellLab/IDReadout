function prototype_readout_forward_model
%% prototype_readout_forward_model.m
% Standalone prototype for MT readout forward model.
%
% Purpose:
%   - Fix an MT population response kernel
%   - Define latent readout functions w(phi)
%   - Predict behavioral scale vs probe offset
%   - Fit Gaussian / Gaussian+offset readout models to observed scales
%
% Current convention:
%   - 0 deg is the task-relevant / preferred direction
%   - Empirical scale is normalized so S_pred(0) = 1 automatically
%   - Fit initially uses only non-anchor offsets (e.g. 45 and 180)
%
% -------------------------------------------------------------------------

clear; clc;

%% ---------------- User inputs ----------------

% Measured pooled scales from kernel analysis
obs.offsetsDeg = [45 180];
obs.scales     = [0.206 0.168];

% Example bootstrap variances; replace with real values when available.
% These are only placeholders for the prototype.
obs.vars       = [0.03^2 0.05^2];

% Fixed MT assumptions
mtParams.sigmaDeg      = 37.5;   % MT direction tuning sigma
mtParams.nullRatioAbs  = 1/3;    % |null response| relative to preferred
                                 % Convention below enforces g(180) = -1/3
mtParams.phiDeg        = -180:1:179;

% Plot settings
plotOffsetsDeg = 0:1:180;

%% ---------------- Build MT kernel ----------------

mt = makeMTKernel(mtParams);

fprintf('MT kernel constants:\n');
fprintf('  sigmaMT = %.2f deg\n', mt.sigmaDeg);
fprintf('  A       = %.6f\n', mt.A);
fprintf('  B       = %.6f\n', mt.B);
fprintf('  g(0)    = %.6f\n', evaluateMTKernel(0, mt));
fprintf('  g(180)  = %.6f\n', evaluateMTKernel(180, mt));
fprintf('\n');

%% ---------------- Fit Gaussian readout ----------------

fitGauss = fitReadoutModel(obs, mt, 'gaussian');

fprintf('Gaussian readout fit:\n');
fprintf('  sigmaW  = %.3f deg\n', fitGauss.params(1));
fprintf('  SSE_wt  = %.6f\n', fitGauss.sseWeighted);
fprintf('  pred(45)= %.6f\n', fitGauss.predAtObs(1));
fprintf('  pred(180)=%.6f\n', fitGauss.predAtObs(2));
fprintf('\n');

%% ---------------- Fit Gaussian + offset readout ----------------

fitGaussOffset = fitReadoutModel(obs, mt, 'gaussian_offset');

fprintf('Gaussian+offset readout fit:\n');
fprintf('  sigmaW  = %.3f deg\n', fitGaussOffset.params(1));
fprintf('  offset  = %.6f\n', fitGaussOffset.params(2));
fprintf('  SSE_wt  = %.6f\n', fitGaussOffset.sseWeighted);
fprintf('  pred(45)= %.6f\n', fitGaussOffset.predAtObs(1));
fprintf('  pred(180)=%.6f\n', fitGaussOffset.predAtObs(2));
fprintf('\n');

%% ---------------- Predict smooth scale curves ----------------

predGauss       = predictScalesFromReadout(plotOffsetsDeg, mt, 'gaussian',        fitGauss.params);
predGaussOffset = predictScalesFromReadout(plotOffsetsDeg, mt, 'gaussian_offset', fitGaussOffset.params);

%% ---------------- Plot 1: MT kernel ----------------

figure(1); clf;
phiPlot = -180:1:180;
plot(phiPlot, evaluateMTKernel(phiPlot, mt), 'k-', 'LineWidth', 1.5);
xlabel('\theta relative to preferred direction (deg)');
ylabel('MT response change, g(\theta)');
title('Fixed MT direction-response kernel');
grid on; box off;
xlim([-180 180]);

%% ---------------- Plot 2: Example readout functions ----------------

figure(2); clf; hold on;
plot(mt.phiDeg, makeReadout(mt.phiDeg, 'gaussian',        fitGauss.params),       'LineWidth', 1.5);
plot(mt.phiDeg, makeReadout(mt.phiDeg, 'gaussian_offset', fitGaussOffset.params), 'LineWidth', 1.5);
xlabel('\phi preferred direction (deg)');
ylabel('Latent readout weight, w(\phi)');
title('Inferred latent readout functions');
legend({'Gaussian','Gaussian + offset'}, 'Location', 'Best');
grid on; box off;
xlim([-180 180]);

%% ---------------- Plot 3: Predicted scale curves and data ----------------

figure(3); clf; hold on;
plot(plotOffsetsDeg, predGauss,       'LineWidth', 1.5);
plot(plotOffsetsDeg, predGaussOffset, 'LineWidth', 1.5);

errorbar(obs.offsetsDeg, obs.scales, sqrt(obs.vars), 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'CapSize', 8);

plot(0, 1, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 7);

xlabel('Probe offset (deg)');
ylabel('Predicted / observed scale');
title('Behavioral scale vs probe offset');
legend({'Gaussian','Gaussian + offset','Observed','Anchor (0,1)'}, ...
    'Location', 'Best');
grid on; box off;
xlim([0 180]);

%% ---------------- Plot 4: MT population drive for selected probe offsets ----------------
% This helps visualize what the readout is integrating.

probeExamples = [0 45 180];
figure(4); clf; hold on;
for k = 1:numel(probeExamples)
    delta = probeExamples(k);
    plot(mt.phiDeg, shiftedMTKernel(mt.phiDeg, delta, mt), 'LineWidth', 1.5);
end
xlabel('\phi preferred direction (deg)');
ylabel('MT population response change');
title('MT population drive for selected probe offsets');
legend(compose('\\Delta = %d^\\circ', probeExamples), 'Location', 'Best');
grid on; box off;
xlim([-180 180]);
end

%% ========================================================================
%% Local functions
%% ========================================================================

function mt = makeMTKernel(params)
% Construct fixed MT directional response-change kernel:
%   g(theta) = A * exp(-theta^2 / (2*sigma^2)) + B
%
% Enforced constraints:
%   g(0)   = 1
%   g(180) = -nullRatioAbs
%
% where theta is wrapped to [-180,180].

    sigmaDeg = params.sigmaDeg;
    phiDeg   = params.phiDeg(:)';
    rNull    = params.nullRatioAbs;

    e180 = exp(-(180^2) / (2 * sigmaDeg^2));

    % Solve:
    %   A + B = 1
    %   A*e180 + B = -rNull
    %
    % => A*(1 - e180) = 1 + rNull
    A = (1 + rNull) / (1 - e180);
    B = 1 - A;

    mt.sigmaDeg = sigmaDeg;
    mt.phiDeg   = phiDeg;
    mt.A        = A;
    mt.B        = B;
end

function g = evaluateMTKernel(thetaDeg, mt)
% Evaluate MT kernel g(theta) with angular wrapping.

    thetaWrapped = wrapTo180Local(thetaDeg);
    g = mt.A .* exp(-(thetaWrapped.^2) ./ (2 * mt.sigmaDeg^2)) + mt.B;
end

function gShift = shiftedMTKernel(phiDeg, deltaDeg, mt)
% MT population response profile across neurons with preferred directions
% phiDeg for a stimulus perturbation at offset deltaDeg:
%   m(phi; delta) = g(phi - delta)

    d = wrapTo180Local(phiDeg - deltaDeg);
    gShift = mt.A .* exp(-(d.^2) ./ (2 * mt.sigmaDeg^2)) + mt.B;
end

function w = makeReadout(phiDeg, modelName, params)
% Construct latent readout function w(phi).
%
% Supported:
%   'gaussian'        params = [sigmaW]
%   'gaussian_offset' params = [sigmaW, offset]
%
% Overall amplitude is irrelevant because prediction is normalized by R(0).

    phiDeg = wrapTo180Local(phiDeg);

    switch lower(modelName)
        case 'gaussian'
            sigmaW = params(1);
            w = exp(-(phiDeg.^2) ./ (2 * sigmaW^2));

        case 'gaussian_offset'
            sigmaW = params(1);
            c      = params(2);
            w = exp(-(phiDeg.^2) ./ (2 * sigmaW^2)) + c;

        otherwise
            error('Unknown modelName: %s', modelName);
    end
end

function pred = predictScalesFromReadout(offsetsDeg, mt, modelName, params)
% Predict normalized behavioral scales:
%   R(delta) = sum_phi w(phi) * m(phi; delta)
%   S_pred(delta) = R(delta) / R(0)

    phiDeg = mt.phiDeg;
    w      = makeReadout(phiDeg, modelName, params);

    R0 = sum(w .* shiftedMTKernel(phiDeg, 0, mt));

    pred = nan(size(offsetsDeg));
    for i = 1:numel(offsetsDeg)
        delta = offsetsDeg(i);
        R = sum(w .* shiftedMTKernel(phiDeg, delta, mt));
        pred(i) = R / R0;
    end
end

function fit = fitReadoutModel(obs, mt, modelName)
% Fit readout parameters to observed non-anchor scales.
%
% Weighted SSE:
%   sum( (obs - pred).^2 ./ var )

    offsetsDeg = obs.offsetsDeg(:)';
    scalesObs  = obs.scales(:)';
    varsObs    = obs.vars(:)';

    switch lower(modelName)
        case 'gaussian'
            % params = [sigmaW]
            p0 = 30;
            lb = 1;
            ub = 180;

        case 'gaussian_offset'
            % params = [sigmaW, offset]
            p0 = [30, 0];
            lb = [1, -2];
            ub = [180, 2];

        otherwise
            error('Unknown modelName: %s', modelName);
    end

    obj = @(p) weightedScaleError(p, offsetsDeg, scalesObs, varsObs, mt, modelName);

    opts = optimset('Display', 'off');
    [pHat, fval] = fminsearchbnd(obj, p0, lb, ub, opts);

    predAtObs = predictScalesFromReadout(offsetsDeg, mt, modelName, pHat);

    fit.modelName    = modelName;
    fit.params       = pHat;
    fit.predAtObs    = predAtObs;
    fit.sseWeighted  = fval;
    fit.offsetsDeg   = offsetsDeg;
    fit.scalesObs    = scalesObs;
    fit.varsObs      = varsObs;
end

function err = weightedScaleError(params, offsetsDeg, scalesObs, varsObs, mt, modelName)

    pred = predictScalesFromReadout(offsetsDeg, mt, modelName, params);
    resid = scalesObs - pred;
    err = sum((resid.^2) ./ varsObs);

    % Keep pathological parameter choices from becoming attractive.
    if ~isfinite(err)
        err = 1e12;
    end
end

function x = wrapTo180Local(x)
% Wrap angle(s) in degrees to [-180, 180].

    x = mod(x + 180, 360) - 180;
end

function [x,fval] = fminsearchbnd(fun,x0,LB,UB,options)
% Simple bounded wrapper around fminsearch.
% This avoids requiring Optimization Toolbox.
%
% Adapted lightweight transform approach for small prototypes.

    x0u = xtransform(x0, LB, UB);
    obj = @(xu) fun(invtransform(xu, LB, UB));
    [xu, fval] = fminsearch(obj, x0u, options);
    x = invtransform(xu, LB, UB);
end

function xu = xtransform(x, LB, UB)
    x = x(:);
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

function x = invtransform(xu, LB, UB)
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