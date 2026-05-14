function toyReadoutDirectionModel
% toyReadoutDirectionModel
%
% Simple toy model for MT tuning + readout weighting.
% Makes:
%   1) Heat map of K(45)/K(0) over readout width and negative offset
%   2) Example weight functions
%   3) Corresponding effective directional sensitivity curves
%
% All angles in degrees.

%% Parameters
sigmaT = 37.5;                    % MT tuning width (deg)
thetaPref = -180:1:180;           % preferred directions for population
thetaStim = -180:1:180;           % stimulus directions for effective kernel

sigmaWVals = 5:1:100;             % readout width sweep
cVals      = 0:0.01:1.00;         % negative offset sweep

% Contour levels for ratio K(45)/K(0)
contourLevels = [-1 -0.5 -0.25 0 0.1 0.20 0.30 0.40 0.50 0.60 0.70 0.75 0.80 0.85 0.90];

%% Precompute MT tuning matrix
% R(s,p) = response of neuron with pref p to stimulus s
R = zeros(numel(thetaStim), numel(thetaPref));
for i = 1:numel(thetaStim)
    d = angdiffDeg(thetaStim(i), thetaPref);
    R(i,:) = exp(-(d.^2) / (2*sigmaT^2));
end

%% Sweep parameter space
ratio45 = nan(numel(cVals), numel(sigmaWVals));
ratio90 = nan(numel(cVals), numel(sigmaWVals)); %#ok<NASGU>
K0mat   = nan(numel(cVals), numel(sigmaWVals));

i0  = find(thetaStim == 0, 1);
i45 = find(thetaStim == 45, 1);
i90 = find(thetaStim == 90, 1);

for ic = 1:numel(cVals)
    c = cVals(ic);
    for isw = 1:numel(sigmaWVals)
      sigmaW = sigmaWVals(isw);
      d0 = angdiffDeg(thetaPref, 0);
      w = exp(-(d0.^2) / (2*sigmaW^2)) - c;

      K = R * w(:);

      K0mat(ic,isw) = K(i0);
      ratio45(ic,isw) = K(i45) / K(i0);
      % ratio90(ic,isw) = K(i90) / K(i0);

      % Signal profile for estimating 0° coherence
      g = exp(-(d0.^2) / (2*sigmaT^2));
      
      % Simple noise model: independent, variance proportional to gain
      var_i = g;   % or ones(size(g)) for equal variance
      signalVal = sum(w .* g);
      noiseVar  = sum((w.^2) .* var_i);
      snrVal = signalVal / sqrt(noiseVar);
      snrMap(ic,isw) = snrVal; %#ok<AGROW>
    end
end

% Mask ill-conditioned points
K0max = max(abs(K0mat(:)), [], 'omitnan');
bad = abs(K0mat) < 0.02 * K0max;   % tune 0.02 as needed
ratio45(bad) = NaN;
% ratio90(bad) = NaN;
% fprintf('min K0 = %.3g, max K0 = %.3g\n', min(K0mat(:)), max(K0mat(:)));
% fprintf('min ratio45 = %.3g, max ratio45 = %.3g\n', min(ratio45(:),[],'omitnan'), max(ratio45(:),[],'omitnan'));
snrMap = snrMap / max(snrMap(:), [], 'omitnan');   % normalize SNR to max = 1

%% Figure 1: heat map + contours for kernel ratio
figure(1); clf
subplot(2,1,1)
imagesc(sigmaWVals, cVals, ratio45, 'AlphaData', ~isnan(ratio45));
axis xy
xlabel('\sigma_W (deg)')
ylabel('negative offset c')
title(sprintf('K(45^\\circ)/K(0^\\circ),  \\sigma_T = %.1f^\\circ', sigmaT))
colorbar
clim([-1 1])   % adjust to taste
hold on
[C,h] = contour(sigmaWVals, cVals, ratio45, contourLevels, 'k', 'LineWidth', 1);
clabel(C,h,'Color','k','FontSize', 9)
[C,h] = contour(sigmaWVals, cVals, ratio45, [0.25 0.25], 'w', 'LineWidth', 1);
clabel(C, h, 'Color', 'w','FontSize', 9);
hold off

%% Figure 1: heat map + contours for SNR
% figure(4); clf
subplot(2,1,2)
imagesc(sigmaWVals, cVals, snrMap);
axis xy
xlabel('\sigma_W (deg)')
ylabel('negative offset c')
title('Relative SNR for 0° coherence estimate')
colorbar
hold on
[C,h] = contour(sigmaWVals, cVals, snrMap, [0.1:0.1:0.9, 0.95, 0.99], 'k', 'LineWidth', 1);
clabel(C, h, 'Color', 'k','FontSize', 9);
[C,h] = contour(sigmaWVals, cVals, ratio45, [0.25 0.25], 'w', 'LineWidth', 1);
clabel(C, h, 'Color', 'w','FontSize', 9);
hold off

% figure(5); clf;
% subplot(1,2,2)
% imagesc(sigmaWVals, cVals, snrMap);
% axis xy
% colorbar
% hold on
% contour(sigmaWVals, cVals, ratio45, [0.25 0.25], 'w', 'LineWidth', 2);
% hold off

%% Figure 2: example weight functions
figure(2); clf

exampleParams = [...
    20 0.00
    40 0.00
    40 0.20
    60 0.35];

nEx = size(exampleParams,1);

for k = 1:nEx
    sigmaW = exampleParams(k,1);
    c      = exampleParams(k,2);

    d0 = angdiffDeg(thetaPref, 0);
    w = exp(-(d0.^2) / (2*sigmaW^2)) - c;

    subplot(2,2,k)
    plot(thetaPref, w, 'LineWidth', 1.5)
    yline(0, ':')
    xlabel('preferred direction (deg)')
    ylabel('weight')
    title(sprintf('\\sigma_W = %g,  c = %.2f', sigmaW, c))
    xlim([-180 180])
end

sgtitle('Example readout weight functions')

%% Figure 3: corresponding effective kernel profiles
figure(3); clf

for k = 1:nEx
    sigmaW = exampleParams(k,1);
    c      = exampleParams(k,2);

    d0 = angdiffDeg(thetaPref, 0);
    w = exp(-(d0.^2) / (2*sigmaW^2)) - c;

    K = R * w(:);
    K = K / K(i0);   % normalize by 0° response

    subplot(2,2,k)
    plot(thetaStim, K, 'LineWidth', 1.5)
    hold on
    plot(0, 1, 'ko', 'MarkerFaceColor', 'k')
    plot(45, K(i45), 'ro', 'MarkerFaceColor', 'r')
    plot(90, K(i90), 'bo', 'MarkerFaceColor', 'b')
    yline(0, ':')
    hold off
    xlabel('stimulus direction (deg)')
    ylabel('K(\theta) / K(0)')
    title(sprintf('\\sigma_W = %g,  c = %.2f', sigmaW, c))
    xlim([-180 180])
end
sgtitle('Effective directional sensitivity')

end

function d = angdiffDeg(a, b)
% Wrapped angular difference a-b in degrees, returned in [-180, 180]
d = mod(a - b + 180, 360) - 180;
end