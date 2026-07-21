%% -------------------------------------------------------------------------
function plotNoiseGainSummary(ax, noiseGain)

hold(ax, 'on');
gainFits = [noiseGain.rect, noiseGain.looKernel];
x = [1 2];
labels = {'Rect', 'LOO'};
colorSpec = {'r', 'b'};
lineSpec = {'-', 'o'};

for i = 1:numel(gainFits)
  gf = gainFits(i);
  plot(ax, [x(i) x(i)], gf.CI95, [colorSpec{i}, lineSpec{1}], 'LineWidth', 1.4, 'HandleVisibility', 'off');
  plot(ax, x(i), gf.gain, [colorSpec{i}, lineSpec{2}], 'MarkerFaceColor', colorSpec{i}, 'MarkerSize', 6, 'DisplayName', labels{i});
end
yline(ax, noiseGain.rect.flatPrediction, 'k--', 'DisplayName', 'flat prediction');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
set(ax, 'XTick', x, 'XTickLabel', labels);
xlim(ax, [0.5 2.5]);
allY = [0, noiseGain.rect.flatPrediction, noiseGain.rect.gain, noiseGain.rect.CI95, ...
  noiseGain.looKernel.gain, noiseGain.looKernel.CI95];
ylimPad = 0.15;
ylim(ax, [min(allY) - ylimPad, max(allY) + ylimPad]);

ylabel(ax, 'Noise Gain');
title(ax, 'Signed-noise Gain');
grid(ax, 'on');
box(ax, 'off');
legend(ax, 'Location', 'best', 'FontSize', 7);

txt = sprintf(['Primary: rect\n' 'rect gain %.3f\n' 'rect z-flat %.2f\n' 'LOO diagnostic\n' ...
  'LOO gain %.3f\n' 'LOO z-flat %.2f'], noiseGain.rect.gain, noiseGain.rect.zVsFlat, ...
  noiseGain.looKernel.gain, noiseGain.looKernel.zVsFlat);
text(ax, 0.03, 0.97, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontSize', 8);

end