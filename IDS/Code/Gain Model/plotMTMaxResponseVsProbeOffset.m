% plotMTMaxResponseVsProbeOffset
% Noiseless, normalized MT population response to symmetric direction probes.
% Streams: 0 deg has a 7% pedestal + 14% step; +/-offset each have 7%.
% Gaussian tuning uses the shortest angular distance around the circle.
% Normalization: response to 100% preferred-direction coherence is 1.
% The full dot mixture contains 35% coherent dots and 65% random dots.
% Random-dot responses are the circular average of the tuning function.
% No offset-dependent rescaling, noise, or additional normalization is used.
% Three additional curves average the response over centered direction
% windows of FULL width 10, 25, and 45 degrees, then maximize over direction.
% Averaging is circular and continuous (Gaussian integrals evaluated by erf).
% Averaged responses retain the original normalization; peaks are not reset.
% Set includeRandomDots=false to plot only coherent-stream contributions.
% This removes a constant, so the curve's shape and maximizing neuron agree.
% Requires only base MATLAB; results remain in the workspace after running.

sigmaDeg = 37.5;
pedestal = 0.07;
step = 0.14;
includeRandomDots = true;
boxcarWidthDeg = [45, 90, 180];

preferredDirectionDeg = (-180:179)'; % 360 distinct preferred directions
probeOffsetDeg = 0:1:180;
coherence = [pedestal + step; pedestal; pedestal];
randomFraction = 1 - sum(coherence);
assert(randomFraction >= 0, 'Total coherence must not exceed 1.');

angularDistance = @(a,b) mod(a - b + 180, 360) - 180;
tuning = @(d) exp(-0.5 .* (d ./ sigmaDeg).^2);

% Antiderivative of periodic, shortest-distance Gaussian tuning. Its linear
% period term lets integration windows cross +/-180 degrees seamlessly.
halfIntegral = sigmaDeg * sqrt(pi/2);
edgeErf = erf(180 / (sqrt(2) * sigmaDeg));
periodIntegral = 2 * halfIntegral * edgeErf;
circularIntegral = @(d) floor((d + 180) ./ 360) .* periodIntegral + ...
    halfIntegral .* (erf(angularDistance(d, 0) ./ (sqrt(2)*sigmaDeg)) + edgeErf);
boxcarTuning = @(d,w) (circularIntegral(d + w/2) - ...
    circularIntegral(d - w/2)) ./ w;

% All preferred directions have the same mean response to random dots.
randomDotResponse = mean(tuning(preferredDirectionDeg));
baseline = double(includeRandomDots) * randomFraction * randomDotResponse;

populationResponse = zeros(numel(preferredDirectionDeg), numel(probeOffsetDeg));
maxResponse = zeros(size(probeOffsetDeg));
maxPreferredDirectionDeg = zeros(size(probeOffsetDeg));
boxcarPopulationResponse = zeros(numel(preferredDirectionDeg), ...
    numel(probeOffsetDeg), numel(boxcarWidthDeg));
maxResponseBoxcar = zeros(numel(boxcarWidthDeg), numel(probeOffsetDeg));
maxPreferredDirectionBoxcarDeg = zeros(size(maxResponseBoxcar));
for k = 1:numel(probeOffsetDeg)
    offset = probeOffsetDeg(k);
    streamDirectionDeg = [0, offset, -offset];
    % bsxfun keeps this compatible with MATLAB versions before R2016b.
    distanceDeg = bsxfun(angularDistance, preferredDirectionDeg, streamDirectionDeg);
    populationResponse(:,k) = baseline + tuning(distanceDeg) * coherence;
    [maxResponse(k), index] = max(populationResponse(:,k));
    maxPreferredDirectionDeg(k) = preferredDirectionDeg(index);
    for j = 1:numel(boxcarWidthDeg)
        % Convolution is linear: average each tuning contribution and sum.
        % The direction-independent random-dot contribution is unchanged.
        boxcarPopulationResponse(:,k,j) = baseline + ...
            boxcarTuning(distanceDeg, boxcarWidthDeg(j)) * coherence;
        [maxResponseBoxcar(j,k), index] = max(boxcarPopulationResponse(:,k,j));
        maxPreferredDirectionBoxcarDeg(j,k) = preferredDirectionDeg(index);
    end
end

figure('Color', 'w', 'Name', 'MT maximum response versus probe offset');
plot(probeOffsetDeg, maxResponse, 'LineWidth', 2.5, 'Color', [0.1 0.35 0.7]);
hold on;
lineColors = [0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];
lineStyles = {'--', '-.', ':'};
for j = 1:numel(boxcarWidthDeg)
    plot(probeOffsetDeg, maxResponseBoxcar(j,:), 'LineWidth', 2.2, ...
        'Color', lineColors(j,:), 'LineStyle', lineStyles{j});
end
hold off;
ax = gca;
set(ax, 'XDir', 'reverse', 'XLim', [0 180], 'XTick', 0:30:180, ...
    'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
xlabel('Symmetric probe offset (degrees)');
ylabel('Maximum normalized population response');
title(sprintf('Gaussian direction tuning: \\sigma = %.1f^{\\circ}', sigmaDeg));
grid on;
legendLabels = [{'No direction averaging'}, ...
    arrayfun(@(w) sprintf('%g-degree boxcar', w), boxcarWidthDeg, ...
    'UniformOutput', false)];
legend(legendLabels, 'Location', 'northwest');

fprintf('Normalization: 100%% preferred-direction coherence = 1.\n');
fprintf('Random-dot contribution: %.6f\n', baseline);
fprintf('Maximum response at 180 deg offset: %.6f\n', maxResponse(end));
fprintf('Maximum response at   0 deg offset: %.6f\n', maxResponse(1));
for j = 1:numel(boxcarWidthDeg)
    fprintf('%2g-degree boxcar: maximum at 180 deg = %.6f; at 0 deg = %.6f\n', ...
        boxcarWidthDeg(j), maxResponseBoxcar(j,end), maxResponseBoxcar(j,1));
end
