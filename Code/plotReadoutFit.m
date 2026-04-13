function plotReadoutFit(fit, mt, varargin)
% plotReadoutFit
%
% Plot:
%   1) inferred latent readout function
%   2) predicted scale curve with observed points
%
% Example:
%   plotReadoutFit(fit, mt, 'titleStr', 'Across-offset readout fit');

p = inputParser;
p.addParameter('titleStr', '', @(x) ischar(x) || isstring(x));
p.addParameter('plotOffsetsDeg', 0:1:180, @(x) validateattributes(x, {'numeric'}, {'vector','real'}));
p.parse(varargin{:});

titleStr       = char(p.Results.titleStr);
plotOffsetsDeg = p.Results.plotOffsetsDeg(:)';

phiDeg = mt.phiDeg;
w      = makeReadout(phiDeg, fit.modelName, fit.params);
pred   = predictScalesFromReadout(plotOffsetsDeg, mt, fit.modelName, fit.params);

figure; clf;

subplot(1,2,1); hold on;
plot(phiDeg, w, 'LineWidth', 1.5);
xlabel('\phi preferred direction (deg)');
ylabel('Readout weight, w(\phi)');
title(sprintf('Latent readout: %s', fit.modelName), 'Interpreter', 'none');
grid on; box off;
xlim([-180 180]);

subplot(1,2,2); hold on;
plot(plotOffsetsDeg, pred, 'LineWidth', 1.5);
errorbar(fit.offsetsFitDeg, fit.obsScale, sqrt(fit.obsVar), 'ko', ...
    'MarkerFaceColor', 'k', 'LineWidth', 1.2, 'CapSize', 8);
plot(0, 1, 'ks', 'MarkerFaceColor', 'w', 'MarkerSize', 7);
xlabel('Probe offset (deg)');
ylabel('Scale');
title('Predicted and observed scale');
grid on; box off;
xlim([0 180]);

if ~isempty(titleStr)
    sgtitle(titleStr);
end

end