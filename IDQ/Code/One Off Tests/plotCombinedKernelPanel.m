%% -------------------------------------------------------------------------
function plotCombinedKernelPanel(ax, kernel)

hold(ax, 'on');

plot(ax, kernel.timeIndex, zeros(size(kernel.timeIndex)), ':', ...
  'HandleVisibility', 'off');

plot(ax, kernel.timeIndex, kernel.meanDiff, '-', ...
  'LineWidth', 1.2, ...
  'DisplayName', 'hit - miss');

title(ax, sprintf('Changed-side signed-noise kernel, n=%d', kernel.nTrials), ...
  'Interpreter', 'none');

xlabel(ax, 'Frame');
ylabel(ax, 'Signed noise, correct - error');
grid(ax, 'on');
box(ax, 'off');

legend(ax, 'Location', 'best', 'FontSize', 7);

txt = sprintf('hit n=%d\nmiss n=%d', kernel.nCorrect, kernel.nError);
text(ax, 0.02, 0.98, txt, ...
  'Units', 'normalized', ...
  'VerticalAlignment', 'top', ...
  'HorizontalAlignment', 'left', ...
  'FontSize', 8);

end