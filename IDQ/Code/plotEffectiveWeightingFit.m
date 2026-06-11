function plotEffectiveWeightingFit(fit, mtModel, varargin)
% plotEffectiveWeightingFit
%
% Plot:
%   1) fitted effective MT-to-choice weighting function
%   2) predicted psychophysical scale curve with observed points
%
% Example:
%   plotEffectiveWeightingFit(fit, mtModel, ...
%       'titleStr', 'Across-offset effective weighting fit');

p = inputParser;
p.addParameter('titleStr', '', @(x) ischar(x) || isstring(x));
p.addParameter('plotOffsetsDeg', 0:1:180, ...
    @(x) validateattributes(x, {'numeric'}, {'vector','real'}));
p.parse(varargin{:});

titleStr       = char(p.Results.titleStr);
plotOffsetsDeg = p.Results.plotOffsetsDeg(:)';

phiDeg = mtModel.phiDeg;

fixedBaseline = [];
if isstruct(fit) && isfield(fit, 'fixedBaseline')
    fixedBaseline = fit.fixedBaseline;
end

wEff = makeEffectiveWeighting(phiDeg, fit.modelName, fit.params, fixedBaseline);
pred = predictScalesFromEffectiveWeighting(plotOffsetsDeg, mtModel, ...
    fit.modelName, fit.params, fixedBaseline);

figure; clf;

subplot(1,2,1); hold on;
plot(phiDeg, wEff, 'LineWidth', 1.5);
xlabel('\phi preferred direction (deg)');
ylabel('Effective weighting, w_{eff}(\phi)');
title(sprintf('Effective weighting: %s', fit.modelName), ...
    'Interpreter', 'none');
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