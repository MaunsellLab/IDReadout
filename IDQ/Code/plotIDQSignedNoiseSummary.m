function plotIDQSignedNoiseSummary()
% plotIDQSignedNoiseSummary
%
% Across-session test for signed-noise asymmetry in the change-side drift
% and pooled non-drift streams. Uses step-rectangular predictors and leaves
% the existing IDQ summaries unchanged.
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_SignedNoiseSummary.mat
%   Plots/AcrossSessionSummaries/IDQSignedNoiseSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQSignedNoiseSummary:NoSessionAnalysisFiles', ...
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

for iDir = 1:3
  sourceName = sprintf('rectStepChangeNoisePredDir%d', iDir);
  trialTable.(sprintf('noisePredDir%d', iDir)) = trialTable.(sourceName);
end

signedFits = fitIDQSignedNoiseAsymmetry(trialTable, ...
  psychDir.sessionFits, alignedWeibull, targetPerformance);
nonDriftCurve = makeBinnedCurve(signedFits.trialData, 'nonDrift', 9);
driftCurve = makeBinnedCurve(signedFits.trialData, 'drift', 9);

T = signedFits.trialData;
streamSD = std([T.nD T.nPlus T.nMinus], 0, 1, 'omitnan');
meanStep = mean(T.stepCoh, 'omitnan');

fprintf('Step-integrated noise SDs: drift %.3f, +120 %.3f, -120 %.3f; mean %.3f%% coherence\n', ...
  streamSD(1), streamSD(2), streamSD(3), mean(streamSD));

fprintf(['Noisy-trial coherence step: mean %.3f%%, SD %.3f%%, median %.3f%%, range %.3f–%.3f%%; ' ...
  'mean stream SD/mean step = %.3f\n'], ...
  meanStep, std(T.stepCoh, 'omitnan'), median(T.stepCoh, 'omitnan'), ...
  min(T.stepCoh), max(T.stepCoh), mean(streamSD) / meanStep);

signedNoiseSummary = struct();
signedNoiseSummary.createdBy = mfilename;
signedNoiseSummary.createdAt = datetime('now');
signedNoiseSummary.noisePredictorType = 'rectStep';
signedNoiseSummary.side = 'change';
signedNoiseSummary.nSessions = numel(sessionAnalyses);
signedNoiseSummary.nTrials = height(trialTable);
signedNoiseSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
signedNoiseSummary.sessionAnalysisFolder = sessionAnalysisFolder;
signedNoiseSummary.initialBetaWeibull = initialBetaWeibull;
signedNoiseSummary.initialLapse = initialLapse;
signedNoiseSummary.targetPerformance = targetPerformance;
signedNoiseSummary.lapseBounds = lapseBounds;
signedNoiseSummary.sessionFits = psychDir.sessionFits;
signedNoiseSummary.alignedWeibull = alignedWeibull;
signedNoiseSummary.trialTable = trialTable;
signedNoiseSummary.signedFits = signedFits;
signedNoiseSummary.nonDriftCurve = nonDriftCurve;
signedNoiseSummary.driftCurve = driftCurve;

summaryFile = fullfile(summaryFolder, 'IDQ_SignedNoiseSummary.mat');
save(summaryFile, 'signedNoiseSummary', '-v7.3');

fig = plotSignedNoiseFigure(signedNoiseSummary);
pdfFile = fullfile(plotFolder, 'IDQSignedNoiseSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
printFitSummary(signedNoiseSummary);

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
    error('plotIDQSignedNoiseSummary:OldSessionAnalysis', ...
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
function curve = makeBinnedCurve(T, streamType, nBins)

switch streamType
  case 'nonDrift'
    x = [T.nPlus; T.nMinus];
    correct = [T.correct; T.correct];
    pMain = [T.pMain; T.pMain];
    pAsymmetry = [T.pNonDriftAsymmetry; T.pNonDriftAsymmetry];
    titleText = 'Individual Non-Drift Streams';
  case 'drift'
    x = T.nD;
    correct = T.correct;
    pMain = T.pMain;
    pAsymmetry = T.pDriftAsymmetry;
    titleText = 'Drift Stream';
  otherwise
    error('plotIDQSignedNoiseSummary:BadStreamType', ...
      'Unknown stream type: %s', streamType);
end

edges = quantile(x, linspace(0, 1, nBins + 1));
edges(1) = -Inf;
edges(end) = Inf;
bin = discretize(x, edges);

meanNoise = nan(nBins, 1);
pObserved = nan(nBins, 1);
pMainBinned = nan(nBins, 1);
pAsymmetryBinned = nan(nBins, 1);
residualMain = nan(nBins, 1);
residualAsymmetry = nan(nBins, 1);
nObservations = zeros(nBins, 1);
for iBin = 1:nBins
  idx = bin == iBin;
  nObservations(iBin) = sum(idx);
  meanNoise(iBin) = mean(x(idx), 'omitnan');
  pObserved(iBin) = mean(correct(idx), 'omitnan');
  pMainBinned(iBin) = mean(pMain(idx), 'omitnan');
  pAsymmetryBinned(iBin) = mean(pAsymmetry(idx), 'omitnan');
  residualMain(iBin) = mean(correct(idx) - pMain(idx), 'omitnan');
  residualAsymmetry(iBin) = mean(correct(idx) - pAsymmetry(idx), 'omitnan');
end

curve = table(meanNoise, pObserved, pMainBinned, pAsymmetryBinned, ...
  residualMain, residualAsymmetry, nObservations);
curve.Properties.Description = titleText;

end

%% ------------------------------------------------------------------------
function fig = plotSignedNoiseFigure(S)

F = S.signedFits;
fig = figure(520);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Change-Side Signed-Noise Asymmetry (%d sessions)', S.nSessions), ...
  'Interpreter', 'none', 'FontWeight', 'bold');

plotChoiceCurve(nexttile(tl, 1), S.nonDriftCurve, 'Non-Drift Choice Relationship');
plotResidualCurve(nexttile(tl, 2), S.nonDriftCurve, 'Non-Drift Residuals');
plotSignedSlopes(nexttile(tl, 3), F.nonDriftSignedSlopes, ...
  'Non-Drift Signed Gains', [0.35 0.55 0.85]);

plotChoiceCurve(nexttile(tl, 4), S.driftCurve, 'Drift Choice Relationship');
plotResidualCurve(nexttile(tl, 5), S.driftCurve, 'Drift Residuals');
plotSignedSlopes(nexttile(tl, 6), F.driftSignedSlopes, ...
  'Drift Signed Gains', [0.85 0.55 0.35]);

plotJointAsymmetry(nexttile(tl, 7), F.bothAsymmetries);
plotDeltaNLL(nexttile(tl, 8), F);
plotTextSummary(nexttile(tl, 9), S);

end

%% ------------------------------------------------------------------------
function plotChoiceCurve(ax, curve, plotTitle)

hold(ax, 'on');
plot(ax, curve.meanNoise, curve.pObserved, 'ko', ...
  'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
plot(ax, curve.meanNoise, curve.pMainBinned, '-', ...
  'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, 'DisplayName', 'Main');
plot(ax, curve.meanNoise, curve.pAsymmetryBinned, '-', ...
  'Color', [0.75 0.25 0.65], 'LineWidth', 1.3, 'DisplayName', 'Signed model');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'Noise predictor');
ylabel(ax, 'P(correct)');
title(ax, plotTitle);
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotResidualCurve(ax, curve, plotTitle)

hold(ax, 'on');
plot(ax, curve.meanNoise, curve.residualMain, 'o-', ...
  'Color', [0.35 0.35 0.35], 'MarkerFaceColor', [0.35 0.35 0.35], ...
  'LineWidth', 1.1, 'DisplayName', 'Main residual');
plot(ax, curve.meanNoise, curve.residualAsymmetry, 'o-', ...
  'Color', [0.75 0.25 0.65], 'MarkerFaceColor', [0.75 0.25 0.65], ...
  'LineWidth', 1.1, 'DisplayName', 'Signed residual');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'Noise predictor');
ylabel(ax, 'Observed - predicted');
title(ax, plotTitle);
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotSignedSlopes(ax, slopes, plotTitle, plotColor)

x = 1:height(slopes);
y = slopes.gain;
low = y - slopes.CI95Low;
high = slopes.CI95High - y;
bar(ax, x, y, 'BarWidth', 0.55, 'FaceColor', plotColor, 'EdgeColor', 'k');
hold(ax, 'on');
errorbar(ax, x, y, low, high, 'k.', 'LineWidth', 1.3);
yline(ax, 0, 'k:');
set(ax, 'XTick', x, 'XTickLabel', slopes.signName);
xlim(ax, [0.25 height(slopes) + 0.75]);
ylabel(ax, 'Noise gain');
title(ax, plotTitle);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotJointAsymmetry(ax, fit)

indices = 3:4;
x = 1:2;
y = fit.gain(indices);
CI = fit.CI95(indices, :);
bar(ax, x, y, 'BarWidth', 0.55, 'FaceColor', [0.65 0.45 0.75], 'EdgeColor', 'k');
hold(ax, 'on');
for i = 1:2
  plot(ax, [x(i) x(i)], CI(i, :), 'k-', 'LineWidth', 1.3);
end
yline(ax, 0, 'k:');
set(ax, 'XTick', x, 'XTickLabel', fit.predictorNames(indices));
xlim(ax, [0.25 2.75]);
ylabel(ax, 'Asymmetry coefficient');
title(ax, 'Joint Asymmetry Terms');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotDeltaNLL(ax, F)

fits = {F.nonDriftAsymmetry, F.driftAsymmetry, F.bothAsymmetries};
labels = {'Non-Drift', 'Drift', 'Both'};
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
function plotTextSummary(ax, S)

F = S.signedFits;
NS = F.nonDriftSignedSlopes;
DS = F.driftSignedSlopes;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s; side: %s', S.noisePredictorType, S.side)
  'Physical coherence signs retained'
  ''
  'Non-drift asymmetry:'
  sprintf('  h_N %.3f [%.3f, %.3f]', F.nonDriftAsymmetry.gain(3), ...
    F.nonDriftAsymmetry.CI95(3, 1), F.nonDriftAsymmetry.CI95(3, 2))
  sprintf('  g_- %.3f; g_+ %.3f', NS.gain(1), NS.gain(2))
  sprintf('  \\DeltaNLL %.3f; p %.3g', F.nonDriftAsymmetry.deltaNLL, ...
    F.nonDriftAsymmetry.pValueLRT)
  ''
  'Drift asymmetry control:'
  sprintf('  h_D %.3f [%.3f, %.3f]', F.driftAsymmetry.gain(3), ...
    F.driftAsymmetry.CI95(3, 1), F.driftAsymmetry.CI95(3, 2))
  sprintf('  g_- %.3f; g_+ %.3f', DS.gain(1), DS.gain(2))
  sprintf('  \\DeltaNLL %.3f; p %.3g', F.driftAsymmetry.deltaNLL, ...
    F.driftAsymmetry.pValueLRT)
  ''
  'g_{positive}=g+h;  g_{negative}=g-h'
  };
text(ax, -0.10, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.signedFits;
NS = F.nonDriftSignedSlopes;
DS = F.driftSignedSlopes;
fprintf('IDQ change-side signed-noise analysis (%s predictor)\n', S.noisePredictorType);
fprintf('Main: gD %.6f, gN %.6f, NLL %.6f\n', ...
  F.main.gain(1), F.main.gain(2), F.main.negLogLikelihood);
fprintf(['Non-drift asymmetry: hN %.6f, gNegative %.6f, gPositive %.6f, ' ...
  'deltaNLL %.6f, LRT p %.6g\n'], ...
  F.nonDriftAsymmetry.gain(3), NS.gain(1), NS.gain(2), ...
  F.nonDriftAsymmetry.deltaNLL, F.nonDriftAsymmetry.pValueLRT);
fprintf(['Drift asymmetry: hD %.6f, gNegative %.6f, gPositive %.6f, ' ...
  'deltaNLL %.6f, LRT p %.6g\n'], ...
  F.driftAsymmetry.gain(3), DS.gain(1), DS.gain(2), ...
  F.driftAsymmetry.deltaNLL, F.driftAsymmetry.pValueLRT);
fprintf('Both: hN %.6f, hD %.6f, deltaNLL %.6f, LRT p %.6g\n', ...
  F.bothAsymmetries.gain(3), F.bothAsymmetries.gain(4), ...
  F.bothAsymmetries.deltaNLL, F.bothAsymmetries.pValueLRT);

end
