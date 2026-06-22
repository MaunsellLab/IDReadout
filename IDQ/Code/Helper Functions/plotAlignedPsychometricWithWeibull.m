%% -------------------------------------------------------------------------
function plotAlignedPsychometricWithWeibull(ax, acrossSummary)

hold(ax, 'on');

T = acrossSummary.trialTable;
fit = acrossSummary.psychFit.pass2.alignedWeibull;

x = T.alignedCoh_pass2;
correct = T.correct;

% Bin aligned coherence for plotting.
binEdges = linspace(min(x), max(x), 9);
binCenters = 0.5 * (binEdges(1:end-1) + binEdges(2:end));

pCorrect = nan(numel(binCenters), 1);
nTrials = nan(numel(binCenters), 1);

for iBin = 1:numel(binCenters)
  idx = x >= binEdges(iBin) & x < binEdges(iBin + 1);

  % Include right edge in final bin.
  if iBin == numel(binCenters)
    idx = x >= binEdges(iBin) & x <= binEdges(iBin + 1);
  end

  nTrials(iBin) = sum(idx);
  if nTrials(iBin) > 0
    pCorrect(iBin) = mean(correct(idx), 'omitnan');
  end
end

idxPlot = nTrials > 0;

plot(ax, binCenters(idxPlot), pCorrect(idxPlot), 'ko', ...
  'MarkerFaceColor', 'k', ...
  'MarkerSize', 4, ...
  'DisplayName', 'binned data');

for i = find(idxPlot(:))'
  text(ax, binCenters(i), pCorrect(i), sprintf(' %d', nTrials(i)), ...
    'FontSize', 7, ...
    'VerticalAlignment', 'bottom');
end

xGrid = linspace(0, max(x) * 1.05, 300);
pFit = idqWeibullP(xGrid, fit.alpha, fit.betaWeibull, fit.lapse);

plot(ax, xGrid, pFit, '-', ...
  'LineWidth', 1.5, ...
  'DisplayName', sprintf('Weibull \\beta=%.2f, lapse=%.3f', ...
  fit.betaWeibull, fit.lapse));

xline(ax, fit.threshold, ':', ...
  'DisplayName', sprintf('threshold%.0f = %.2f', ...
  100 * fit.thresholdPerformance, fit.threshold));

yline(ax, fit.thresholdPerformance, ':', ...
  'HandleVisibility', 'off');

xlabel(ax, sprintf('Coherence / session threshold%.0f', ...
  100 * fit.thresholdPerformance));
ylabel(ax, 'P(correct)');
title(ax, 'Aligned psychometric with Weibull fit');

ylim(ax, [0.45 1.02]);
xlim(ax, [0 max(xGrid)]);
grid(ax, 'on');
box(ax, 'off');
legend(ax, 'Location', 'best', 'FontSize', 7);

end