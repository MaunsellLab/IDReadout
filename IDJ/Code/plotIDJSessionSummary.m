function [summary, fig] = plotIDJSessionSummary(source, varargin)
% plotIDJSessionSummary  Plot one compact IDJ session on a single page.
%
% source can be a sessionAnalysis structure or a *_sessionAnalysis.mat file.

p = inputParser;
p.addParameter('LapseBounds', [0 0.15], @(x) isnumeric(x) && numel(x) == 2);
% PDFPath is retained for compatibility, but only its filename is used.
% Every PDF is written beneath IDJ/Plots/FullSessions.
p.addParameter('PDFPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('ExportPDF', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Visible', 'on', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

if ischar(source) || isstring(source)
  loaded = load(source, 'sessionAnalysis');
  A = loaded.sessionAnalysis;
else
  A = source;
end

assert(idjHasTopUpExclusion(A), 'plotIDJSessionSummary:RebuildAnalysis', ...
  'Rebuild compact sessions with makeIDJSessionAnalyses(''Replace'',true); top-up exclusion is required.');
T = A.trialTable;
freeLapse = fitIDJSessionPsychometrics(T, 'LapseBounds', opts.LapseBounds);
zeroLapse = fitIDJSessionPsychometrics(T, 'FixedLapse', 0, ...
  'LapseBounds', opts.LapseBounds);

summary = struct();
summary.sourceSessionID = A.sourceSessionID;
summary.createdBy = mfilename;
summary.createdAt = datetime('now');
summary.freeLapseFit = freeLapse;
summary.zeroLapseFit = zeroLapse;
summary.configurationTable = A.configurationTable;
summary.stepTable = A.stepTable;
summary.pdfPath = "";

fig = figure('Color', 'w', 'Visible', char(opts.Visible), ...
  'Position', [100 100 1000 900]);
tl = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('IDJ Session Summary: %s', char(A.sourceSessionID)), ...
  'Interpreter', 'none', 'FontWeight', 'bold');

allocationPC = freeLapse.thresholdTable.allocationPC;
colors = [0.10 0.45 0.80; 0.90 0.45 0.10; 0.55 0.30 0.70; 0.20 0.65 0.40];
markers = {'o','s','^','d','v','>'};
maxStep = max(T.totalStepCohPC);
for k = 1:numel(allocationPC)
  ax = nexttile(tl, k);
  plotOneAllocation(ax, T, zeroLapse, allocationPC(k), colors(k,:), ...
    colors, markers, maxStep);
end

plotThresholdSummary(nexttile(tl, 5), freeLapse, zeroLapse, colors);
plotSessionText(nexttile(tl, 6), A, freeLapse, zeroLapse);

if opts.ExportPDF
  if strlength(string(opts.PDFPath)) == 0
    pdfName = sprintf('%s_IDJSummary.pdf', char(A.sourceSessionID));
  else
    [~, pdfBase, pdfExtension] = fileparts(char(opts.PDFPath));
    if isempty(pdfExtension)
      pdfExtension = '.pdf';
    end
    assert(strcmpi(pdfExtension, '.pdf'), ...
      'plotIDJSessionSummary:PDFExtension', ...
      'The requested graphical output must be a PDF.');
    pdfName = [pdfBase pdfExtension];
  end
  pdfPath = fullfile(idjPlotFolder('FullSessions'), pdfName);
  exportgraphics(fig, pdfPath, 'ContentType', 'vector');
  summary.pdfPath = string(pdfPath);
end
end

%% ------------------------------------------------------------------------
function plotOneAllocation(ax, T, fit, allocationPC, color, colors, markers, maxStep)

hold(ax, 'on');
allocationMask = abs(T.unsignedAllocPC - allocationPC) < 1e-6;
epochs = unique(T.configurationEpoch(allocationMask), 'stable');
for e = 1:numel(epochs)
  epochMask = allocationMask & T.configurationEpoch == epochs(e);
  steps = unique(T.totalStepCohPC(epochMask));
  for s = 1:numel(steps)
    idx = epochMask & T.totalStepCohPC == steps(s);
    proportion = mean(T.correct(idx));
    plot(ax, steps(s), proportion, markers{1 + mod(e - 1, numel(markers))}, ...
      'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 5, ...
      'HandleVisibility', 'off');
    verticalAlignment = 'bottom';
    if proportion > 0.96
      verticalAlignment = 'top';
    end
    text(ax, steps(s), proportion, sprintf(' %d', sum(idx)), ...
      'FontSize', 7, 'VerticalAlignment', verticalAlignment, ...
      'Color', color, 'Clipping', 'on');
  end
end

row = find(abs(fit.thresholdTable.allocationPC - allocationPC) < 1e-6);
x = linspace(0, maxStep * 1.08, 300);
y = idjWeibullP(x, fit.thresholdTable.alpha(row), fit.betaWeibull, fit.lapse);
plot(ax, x, y, '-', 'Color', color, 'LineWidth', 2);

for k = 1:height(fit.thresholdTable)
  width = 2.0;
  if k == row
    width = 2.0;
  end
  xline(ax, fit.thresholdTable.threshold(k), ':', 'Color', colors(k,:), ...
    'LineWidth', width, 'HandleVisibility', 'off');
end
yline(ax, fit.targetPerformance, ':', 'Color', [0.35 0.35 0.35], ...
  'HandleVisibility', 'off');
text(ax, 0.03, 0.96, sprintf('75%% threshold: %.1f%%', ...
  fit.thresholdTable.threshold(row)), 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'FontSize', 8, 'FontWeight', 'bold', ...
  'Color', color, 'BackgroundColor', 'w', ...
  'EdgeColor', [0.80 0.80 0.80], 'Margin', 2, 'Clipping', 'on');
title(ax, allocationLabel(allocationPC));
xlabel(ax, 'Total coherence step (%)');
ylabel(ax, 'P(correct)');
xlim(ax, [0 maxStep * 1.08]);
ylim(ax, [0.45 1.02]);
grid(ax, 'on');
box(ax, 'off');
end

%% ------------------------------------------------------------------------
function plotThresholdSummary(ax, freeFit, zeroFit, colors)

hold(ax, 'on');
x = (1:height(zeroFit.thresholdTable))';
threshold = zeroFit.thresholdTable.threshold;
low = zeroFit.thresholdTable.thresholdCILow;
high = zeroFit.thresholdTable.thresholdCIHigh;
lowerError = threshold - low;
upperError = high - threshold;
for k = 1:numel(x)
  errorbar(ax, x(k), threshold(k), lowerError(k), upperError(k), 'o', ...
    'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), 'LineWidth', 1.1, ...
    'HandleVisibility', 'off');
end
plot(ax, x, freeFit.thresholdTable.threshold, 'ko', ...
  'MarkerFaceColor', 'none', 'DisplayName', 'Free shared lapse');
plot(ax, nan, nan, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'Zero lapse');
set(ax, 'XTick', x, 'XTickLabel', arrayfun(@allocationLabel, zeroFit.thresholdTable.allocationPC, 'UniformOutput', false));
xlabel(ax, 'Allocation');
ylabel(ax, '75% threshold (%)');
title(ax, 'Threshold Comparison');
legend(ax, 'Location', 'best');
grid(ax, 'on');
box(ax, 'off');
end

%% ------------------------------------------------------------------------
function plotSessionText(ax, A, freeFit, zeroFit)

axis(ax, 'off');
H = A.sessionHeader;
lines = { ...
  sprintf('Animal: %s', H.animal); ...
  sprintf('Date: %s', char(string(H.date))); ...
  sprintf('All trials: %d', A.nTrialsTotal); ...
  sprintf('Analysis trials: %d', A.nAnalysisTrials); ...
  sprintf('Excluded top-ups: %d', A.nExcludedTopUpTrials); ...
  sprintf('Correct: %.1f%%', 100 * A.nCorrect / A.nAnalysisTrials); ...
  sprintf('Configuration epochs: %d', height(A.configurationTable)); ...
  ''; ...
  sprintf('Zero-lapse beta: %.3f', zeroFit.betaWeibull); ...
  sprintf('Free-lapse beta: %.3f', freeFit.betaWeibull); ...
  sprintf('Shared lapse: %.3f [%.3f, %.3f]', ...
    freeFit.lapse, freeFit.lapseCI95(1), freeFit.lapseCI95(2)); ...
  sprintf('Free-lapse NLL: %.1f', freeFit.negLogLikelihood); ...
  sprintf('Zero-lapse NLL: %.1f', zeroFit.negLogLikelihood); ...
  sprintf('2 Delta NLL: %.1f', ...
    2 * (zeroFit.negLogLikelihood - freeFit.negLogLikelihood)); ...
  ''; ...
  'Point labels give trial counts.'; ...
  'Point symbols distinguish configuration epochs.' ...
  };
text(ax, 0, 1, lines, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
  'FontName', 'Menlo', 'FontSize', 8);
end

%% ------------------------------------------------------------------------
function label = allocationLabel(majorPC)
if abs(majorPC - 100/3) < 1e-6
  label = '33.3:33.3:33.3';
else
  label = sprintf('%g:%g:%g', majorPC, (100 - majorPC)/2, (100 - majorPC)/2);
end
end
