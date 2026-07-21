function plotIDQInteractionSummary()
% plotIDQInteractionSummary
%
% Focused across-session analysis of interactions among the three
% change-side IDQ noise streams. This function leaves plotIDQSummary and
% plotIDQSideSummary unchanged.
%
% The current analysis uses step-rectangular predictors and fits:
%   main: gD*nD + gN*(nPlus+nMinus)
%   main + hDN*nD*(nPlus+nMinus)
%   main + hNN*nPlus*nMinus
%   main + both interaction terms
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_InteractionSummary.mat
%   Plots/AcrossSessionSummaries/IDQInteractionSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQInteractionSummary:NoSessionAnalysisFiles', ...
    'No sessionAnalysis files found in %s.', sessionAnalysisFolder);
end

sessionAnalyses = cell(numel(files), 1);
for iFile = 1:numel(files)
  filePath = fullfile(files(iFile).folder, files(iFile).name);
  load(filePath, 'sessionAnalysis');
  sessionAnalyses{iFile} = sessionAnalysis;
end
validateInputs(sessionAnalyses);

trialTable = concatenateTrialTables(sessionAnalyses);
initialBetaWeibull = 2;
initialLapse = 0.02;
targetPerformance = 0.75;
lapseBounds = [0 0.05];

psychDir = fitIDQInitialSessionThresholds(trialTable, ...
  initialBetaWeibull, initialLapse, targetPerformance);
alignedWeibull = fitIDQAcrossAlignedWeibull(psychDir.alignedCoh, ...
  trialTable.correct, targetPerformance, lapseBounds);
trialTable.alignedCoh = psychDir.alignedCoh;
trialTable.noisyStepAlignedCoh = psychDir.noisyStepAlignedCoh;

% Use the change-side step-rectangular streams without changing their signs.
for iDir = 1:3
  sourceName = sprintf('rectStepChangeNoisePredDir%d', iDir);
  trialTable.(sprintf('noisePredDir%d', iDir)) = trialTable.(sourceName);
end

interactionFits = fitIDQNoiseInteractions(trialTable, ...
  psychDir.sessionFits, alignedWeibull, targetPerformance);
surfaceSummary = computeBinnedSurfaces(interactionFits.trialData, 4);

interactionSummary = struct();
interactionSummary.createdBy = mfilename;
interactionSummary.createdAt = datetime('now');
interactionSummary.noisePredictorType = 'rectStep';
interactionSummary.side = 'change';
interactionSummary.nSessions = numel(sessionAnalyses);
interactionSummary.nTrials = height(trialTable);
interactionSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
interactionSummary.sessionAnalysisFolder = sessionAnalysisFolder;
interactionSummary.initialBetaWeibull = initialBetaWeibull;
interactionSummary.initialLapse = initialLapse;
interactionSummary.targetPerformance = targetPerformance;
interactionSummary.lapseBounds = lapseBounds;
interactionSummary.sessionFits = psychDir.sessionFits;
interactionSummary.alignedWeibull = alignedWeibull;
interactionSummary.trialTable = trialTable;
interactionSummary.interactionFits = interactionFits;
interactionSummary.surfaceSummary = surfaceSummary;

summaryFile = fullfile(summaryFolder, 'IDQ_InteractionSummary.mat');
save(summaryFile, 'interactionSummary', '-v7.3');

fig = plotInteractionFigure(interactionSummary);
pdfFile = fullfile(plotFolder, 'IDQInteractionSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

printFitSummary(interactionSummary);

end

%% ------------------------------------------------------------------------
function validateInputs(sessionAnalyses)

requiredVars = strings(1, 3);
for iDir = 1:3
  requiredVars(iDir) = sprintf('rectStepChangeNoisePredDir%d', iDir);
end

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  missing = requiredVars(~ismember(requiredVars, ...
    string(SA.trialTable.Properties.VariableNames)));
  if ~isempty(missing)
    error('plotIDQInteractionSummary:OldSessionAnalysis', ...
      ['Session %s lacks the change-side rectangular predictors. Rerun ' ...
       'makeIDQSessionAnalyses before this analysis.'], SA.fileName);
  end
end

end

%% ------------------------------------------------------------------------
function trialTable = concatenateTrialTables(sessionAnalyses)

tables = cell(numel(sessionAnalyses), 1);
for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  T = SA.trialTable;
  T.sessionIndex = repmat(iSession, height(T), 1);
  T.noisyStepCoh = repmat(SA.noisyStepCoh, height(T), 1);
  tables{iSession} = T;
end
trialTable = vertcat(tables{:});
trialTable = movevars(trialTable, 'sessionIndex', 'Before', 1);

end

%% ------------------------------------------------------------------------
function surfaceSummary = computeBinnedSurfaces(T, nBins)

edgesD = quantile(T.nD, linspace(0, 1, nBins + 1));
edgesN = quantile(T.nN, linspace(0, 1, nBins + 1));
edgesD(1) = -Inf;
edgesD(end) = Inf;
edgesN(1) = -Inf;
edgesN(end) = Inf;

binD = discretize(T.nD, edgesD);
binN = discretize(T.nN, edgesN);

observed = nan(nBins);
predictedMain = nan(nBins);
predictedBoth = nan(nBins);
residualMain = nan(nBins);
nTrials = zeros(nBins);
meanND = nan(nBins);
meanNN = nan(nBins);

% Rows are pooled non-drift bins; columns are drift bins.
for iN = 1:nBins
  for iD = 1:nBins
    idx = binD == iD & binN == iN;
    nTrials(iN, iD) = sum(idx);
    if any(idx)
      observed(iN, iD) = mean(T.correct(idx), 'omitnan');
      predictedMain(iN, iD) = mean(T.pMain(idx), 'omitnan');
      predictedBoth(iN, iD) = mean(T.pBoth(idx), 'omitnan');
      residualMain(iN, iD) = mean(T.correct(idx) - T.pMain(idx), 'omitnan');
      meanND(iN, iD) = mean(T.nD(idx), 'omitnan');
      meanNN(iN, iD) = mean(T.nN(idx), 'omitnan');
    end
  end
end

surfaceSummary = struct();
surfaceSummary.nBins = nBins;
surfaceSummary.edgesD = edgesD;
surfaceSummary.edgesN = edgesN;
surfaceSummary.binD = binD;
surfaceSummary.binN = binN;
surfaceSummary.observed = observed;
surfaceSummary.predictedMain = predictedMain;
surfaceSummary.predictedBoth = predictedBoth;
surfaceSummary.residualMain = residualMain;
surfaceSummary.nTrials = nTrials;
surfaceSummary.meanND = meanND;
surfaceSummary.meanNN = meanNN;

end

%% ------------------------------------------------------------------------
function fig = plotInteractionFigure(S)

F = S.interactionFits;
surf = S.surfaceSummary;

fig = figure(510);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Change-Side Noise Interactions (%d sessions)', S.nSessions), ...
  'Interpreter', 'none', 'FontWeight', 'bold');

probValues = [surf.observed(:); surf.predictedMain(:)];
probValues = probValues(isfinite(probValues));
probLim = [min(probValues), max(probValues)];
if diff(probLim) == 0
  probLim = probLim + [-0.01 0.01];
end

plotSurface(nexttile(tl, 1), surf.observed, probLim, ...
  'Observed P(correct)', surf.nTrials, false);
plotSurface(nexttile(tl, 2), surf.predictedMain, probLim, ...
  'Main-Effects Prediction', [], false);

residMax = max(abs(surf.residualMain(:)), [], 'omitnan');
if ~isfinite(residMax) || residMax == 0
  residMax = 0.01;
end
plotSurface(nexttile(tl, 3), surf.residualMain, [-residMax residMax], ...
  'Observed - Main Prediction', [], true);

plotCoefficientBars(nexttile(tl, 4), F.main, 1:2, ...
  'Reduced Main Effects', [0.35 0.55 0.85]);
plotCoefficientBars(nexttile(tl, 5), F.bothInteractions, 3:4, ...
  'Joint Interaction Terms', [0.75 0.45 0.85]);
plotDeltaNLL(nexttile(tl, 6), F);
plotStratifiedGain(nexttile(tl, 7), F.stratifiedMain);
plotConditionalResiduals(nexttile(tl, 8), F.trialData);
plotTextSummary(nexttile(tl, 9), S);

end


%% ------------------------------------------------------------------------
function plotSurface(ax, values, colorLimits, plotTitle, nTrials, diverging)

imagesc(ax, values, colorLimits);
set(ax, 'YDir', 'normal', 'XTick', 1:size(values, 2), ...
  'YTick', 1:size(values, 1));
xlabel(ax, 'Drift-noise quartile');
ylabel(ax, 'Pooled non-drift quartile');
title(ax, plotTitle, 'Interpreter', 'none');
box(ax, 'off');
colorbar(ax);
if diverging
  colormap(ax, blueWhiteRed(256));
else
  colormap(ax, parula(256));
end

for iRow = 1:size(values, 1)
  for iCol = 1:size(values, 2)
    if isfinite(values(iRow, iCol))
      if isempty(nTrials)
        txt = sprintf('%.3f', values(iRow, iCol));
      else
        txt = sprintf('%.3f\nn=%d', values(iRow, iCol), nTrials(iRow, iCol));
      end
      text(ax, iCol, iRow, txt, 'HorizontalAlignment', 'center', ...
        'FontSize', 7, 'Color', contrastTextColor(values(iRow, iCol), colorLimits));
    end
  end
end

end

%% ------------------------------------------------------------------------
function color = contrastTextColor(value, limits)

relative = (value - limits(1)) / max(diff(limits), eps);
if relative < 0.25 || relative > 0.80
  color = 'w';
else
  color = 'k';
end

end

%% ------------------------------------------------------------------------
function cmap = blueWhiteRed(n)

if nargin < 1
  n = 256;
end
nLower = floor(n / 2);
nUpper = n - nLower;
cmap = [linspace(0.1, 1, nLower)', linspace(0.3, 1, nLower)', ones(nLower, 1); ...
  ones(nUpper, 1), linspace(1, 0.2, nUpper)', linspace(1, 0.2, nUpper)'];

end

%% ------------------------------------------------------------------------
function plotCoefficientBars(ax, fit, indices, plotTitle, plotColor)

g = fit.gain(indices);
CI = fit.CI95(indices, :);
x = 1:numel(indices);
bar(ax, x, g, 'BarWidth', 0.55, 'FaceColor', plotColor, 'EdgeColor', 'k');
hold(ax, 'on');
for i = 1:numel(indices)
  plot(ax, [x(i) x(i)], CI(i, :), 'k-', 'LineWidth', 1.3);
end
plot(ax, x, g, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
  'LineStyle', 'none');
yline(ax, 0, 'k:');
set(ax, 'XTick', x, 'XTickLabel', fit.predictorNames(indices));
xlim(ax, [0.25, numel(indices) + 0.75]);
ylabel(ax, 'Coefficient');
title(ax, plotTitle, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotDeltaNLL(ax, F)

fits = {F.driftNonDriftInteraction, F.nonDriftInteraction, F.bothInteractions};
labels = {'h_{DN}', 'h_{NN}', 'Both'};
deltaNLL = cellfun(@(x) x.deltaNLL, fits);
bar(ax, 1:3, deltaNLL, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'k');
set(ax, 'XTick', 1:3, 'XTickLabel', labels);
ylabel(ax, '\DeltaNLL from main model');
title(ax, 'Nested Model Improvement');
grid(ax, 'on');
box(ax, 'off');
for i = 1:3
  text(ax, i, deltaNLL(i), sprintf('p=%.3g', fits{i}.pValueLRT), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 8);
end

end

%% ------------------------------------------------------------------------
function plotStratifiedGain(ax, stratified)

x = stratified.medianND;
y = stratified.gN;
low = y - stratified.gNCI95Low;
high = stratified.gNCI95High - y;
errorbar(ax, x, y, low, high, 'ko-', 'MarkerFaceColor', 'k', ...
  'LineWidth', 1.2, 'MarkerSize', 5);
yline(ax, 0, 'k:');
xlabel(ax, 'Median drift-noise predictor');
ylabel(ax, 'Pooled non-drift gain');
title(ax, 'Non-Drift Gain by Drift Stratum');
grid(ax, 'on');
box(ax, 'off');
for i = 1:height(stratified)
  text(ax, x(i), stratified.gNCI95High(i), sprintf(' n=%d', stratified.nTrials(i)), ...
    'FontSize', 7, 'VerticalAlignment', 'bottom');
end

end

%% ------------------------------------------------------------------------
function plotConditionalResiduals(ax, T)

nStrata = 3;
nBins = 5;
edgesD = quantile(T.nD, [0, 1/3, 2/3, 1]);
edgesD(1) = -Inf;
edgesD(end) = Inf;
colors = lines(nStrata);
hold(ax, 'on');

for iStratum = 1:nStrata
  idxStratum = T.nD >= edgesD(iStratum) & T.nD < edgesD(iStratum + 1);
  if iStratum == nStrata
    idxStratum = T.nD >= edgesD(iStratum) & T.nD <= edgesD(iStratum + 1);
  end
  nNThis = T.nN(idxStratum);
  residualThis = T.correct(idxStratum) - T.pMain(idxStratum);
  edgesN = quantile(nNThis, linspace(0, 1, nBins + 1));
  edgesN(1) = -Inf;
  edgesN(end) = Inf;
  binN = discretize(nNThis, edgesN);
  meanN = nan(nBins, 1);
  meanResidual = nan(nBins, 1);
  for iBin = 1:nBins
    idxBin = binN == iBin;
    meanN(iBin) = mean(nNThis(idxBin), 'omitnan');
    meanResidual(iBin) = mean(residualThis(idxBin), 'omitnan');
  end
  plot(ax, meanN, meanResidual, 'o-', 'Color', colors(iStratum, :), ...
    'MarkerFaceColor', colors(iStratum, :), 'LineWidth', 1.1, ...
    'DisplayName', sprintf('Drift stratum %d', iStratum));
end

yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'Pooled non-drift predictor');
ylabel(ax, 'Observed - main P(correct)');
title(ax, 'Conditional Main-Model Residuals');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotTextSummary(ax, S)

F = S.interactionFits;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s; side: %s', S.noisePredictorType, S.side)
  'Physical coherence signs retained'
  ''
  'Main:'
  sprintf('  g_D %.3f [%.3f, %.3f]', F.main.gain(1), F.main.CI95(1, 1), F.main.CI95(1, 2))
  sprintf('  g_N %.3f [%.3f, %.3f]', F.main.gain(2), F.main.CI95(2, 1), F.main.CI95(2, 2))
  sprintf('  NLL %.3f', F.main.negLogLikelihood)
  ''
  'Joint interactions:'
  sprintf('  h_DN %.3f [%.3f, %.3f]', F.bothInteractions.gain(3), ...
    F.bothInteractions.CI95(3, 1), F.bothInteractions.CI95(3, 2))
  sprintf('  h_NN %.3f [%.3f, %.3f]', F.bothInteractions.gain(4), ...
    F.bothInteractions.CI95(4, 1), F.bothInteractions.CI95(4, 2))
  sprintf('  \DeltaNLL %.3f; LRT p %.3g', F.bothInteractions.deltaNLL, ...
    F.bothInteractions.pValueLRT)
  ''
  'c_{eff}=c_{step}+g_Dn_D+g_N(n_++n_-)'
  '       +h_{DN}I_{DN}+h_{NN}I_{NN}'
  };
text(ax, -0.10, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.interactionFits;
fprintf('IDQ change-side interaction analysis (%s predictor)\n', S.noisePredictorType);
fprintf('Main: gD %.6f, gN %.6f, NLL %.6f\n', ...
  F.main.gain(1), F.main.gain(2), F.main.negLogLikelihood);
fprintf('Add hDN: hDN %.6f, deltaNLL %.6f, LRT p %.6g\n', ...
  F.driftNonDriftInteraction.gain(3), ...
  F.driftNonDriftInteraction.deltaNLL, ...
  F.driftNonDriftInteraction.pValueLRT);
fprintf('Add hNN: hNN %.6f, deltaNLL %.6f, LRT p %.6g\n', ...
  F.nonDriftInteraction.gain(3), F.nonDriftInteraction.deltaNLL, ...
  F.nonDriftInteraction.pValueLRT);
fprintf('Both: hDN %.6f, hNN %.6f, deltaNLL %.6f, LRT p %.6g\n', ...
  F.bothInteractions.gain(3), F.bothInteractions.gain(4), ...
  F.bothInteractions.deltaNLL, F.bothInteractions.pValueLRT);

end
