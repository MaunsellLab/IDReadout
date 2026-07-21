function plotIDQNoChangeSignedNoiseSummary()
% plotIDQNoChangeSignedNoiseSummary
%
% Test signed-noise asymmetry across the three no-change-side streams using
% a common gain and a common asymmetry term. Uses step-rectangular predictors.
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_NoChangeSignedNoiseSummary.mat
%   Plots/AcrossSessionSummaries/IDQNoChangeSignedNoiseSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQNoChangeSignedNoiseSummary:NoSessionAnalysisFiles', ...
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
  sourceName = sprintf('rectStepNoChangeNoisePredDir%d', iDir);
  trialTable.(sprintf('noisePredDir%d', iDir)) = trialTable.(sourceName);
end

signedFits = fitIDQNoChangeSignedNoiseAsymmetry(trialTable, ...
  psychDir.sessionFits, alignedWeibull, targetPerformance);
allStreamCurve = makeBinnedCurve(signedFits.trialData, 'all', 9);
relativeCurves = cell(1, 3);
relativeNames = {'Drift-matched', '+120°', '-120°'};
for iDir = 1:3
  relativeCurves{iDir} = makeBinnedCurve(signedFits.trialData, iDir, 9);
end

noChangeSignedNoiseSummary = struct();
noChangeSignedNoiseSummary.createdBy = mfilename;
noChangeSignedNoiseSummary.createdAt = datetime('now');
noChangeSignedNoiseSummary.noisePredictorType = 'rectStep';
noChangeSignedNoiseSummary.side = 'noChange';
noChangeSignedNoiseSummary.nSessions = numel(sessionAnalyses);
noChangeSignedNoiseSummary.nTrials = height(trialTable);
noChangeSignedNoiseSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
noChangeSignedNoiseSummary.sessionAnalysisFolder = sessionAnalysisFolder;
noChangeSignedNoiseSummary.initialBetaWeibull = initialBetaWeibull;
noChangeSignedNoiseSummary.initialLapse = initialLapse;
noChangeSignedNoiseSummary.targetPerformance = targetPerformance;
noChangeSignedNoiseSummary.lapseBounds = lapseBounds;
noChangeSignedNoiseSummary.sessionFits = psychDir.sessionFits;
noChangeSignedNoiseSummary.alignedWeibull = alignedWeibull;
noChangeSignedNoiseSummary.trialTable = trialTable;
noChangeSignedNoiseSummary.signedFits = signedFits;
noChangeSignedNoiseSummary.allStreamCurve = allStreamCurve;
noChangeSignedNoiseSummary.relativeCurves = relativeCurves;
noChangeSignedNoiseSummary.relativeNames = relativeNames;

summaryFile = fullfile(summaryFolder, 'IDQ_NoChangeSignedNoiseSummary.mat');
save(summaryFile, 'noChangeSignedNoiseSummary', '-v7.3');

fig = plotSignedNoiseFigure(noChangeSignedNoiseSummary);
pdfFile = fullfile(plotFolder, 'IDQNoChangeSignedNoiseSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
printFitSummary(noChangeSignedNoiseSummary);

end

%% ------------------------------------------------------------------------
function validateInputs(sessionAnalyses)

requiredVars = strings(1, 3);
for iDir = 1:3
  requiredVars(iDir) = sprintf('rectStepNoChangeNoisePredDir%d', iDir);
end
for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  missing = requiredVars(~ismember(requiredVars, ...
    string(SA.trialTable.Properties.VariableNames)));
  if ~isempty(missing)
    error('plotIDQNoChangeSignedNoiseSummary:OldSessionAnalysis', ...
      ['Session %s lacks the no-change-side rectangular predictors. Rerun ' ...
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
function curve = makeBinnedCurve(T, streamSelection, nBins)

if isequal(streamSelection, 'all')
  x = [T.nD; T.nPlus; T.nMinus];
  correct = repmat(T.correct, 3, 1);
  pMain = repmat(T.pMain, 3, 1);
  pSigned = repmat(T.pSigned, 3, 1);
else
  streamNames = {'nD', 'nPlus', 'nMinus'};
  x = T.(streamNames{streamSelection});
  correct = T.correct;
  pMain = T.pMain;
  pSigned = T.pSigned;
end

edges = quantile(x, linspace(0, 1, nBins + 1));
edges(1) = -Inf;
edges(end) = Inf;
bin = discretize(x, edges);
meanNoise = nan(nBins, 1);
pObserved = nan(nBins, 1);
pMainBinned = nan(nBins, 1);
pSignedBinned = nan(nBins, 1);
pOpponentBinned = nan(nBins, 1);
residualMain = nan(nBins, 1);
residualSigned = nan(nBins, 1);
residualOpponent = nan(nBins, 1);
nObservations = zeros(nBins, 1);
for iBin = 1:nBins
  idx = bin == iBin;
  nObservations(iBin) = sum(idx);
  meanNoise(iBin) = mean(x(idx), 'omitnan');
  pObserved(iBin) = mean(correct(idx), 'omitnan');
  pMainBinned(iBin) = mean(pMain(idx), 'omitnan');
  pSignedBinned(iBin) = mean(pSigned(idx), 'omitnan');
  if isequal(streamSelection, 'all')
    pOpponent = repmat(T.pOpponentConstrained, 3, 1);
  else
    pOpponent = T.pOpponentConstrained;
  end
  pOpponentBinned(iBin) = mean(pOpponent(idx), 'omitnan');
  residualMain(iBin) = mean(correct(idx) - pMain(idx), 'omitnan');
  residualSigned(iBin) = mean(correct(idx) - pSigned(idx), 'omitnan');
  residualOpponent(iBin) = mean(correct(idx) - pOpponent(idx), 'omitnan');
end
curve = table(meanNoise, pObserved, pMainBinned, pSignedBinned, ...
  pOpponentBinned, residualMain, residualSigned, residualOpponent, ...
  nObservations);

end

%% ------------------------------------------------------------------------
function fig = plotSignedNoiseFigure(S)

F = S.signedFits;
fig = figure(530);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ No-Change-Side Signed-Noise Asymmetry (%d sessions)', S.nSessions), ...
  'Interpreter', 'none', 'FontWeight', 'bold');

plotChoiceCurve(nexttile(tl, 1), S.allStreamCurve);
plotResidualCurve(nexttile(tl, 2), S.allStreamCurve);
plotSignedSlopes(nexttile(tl, 3), F.signedSlopes, F.opponentConstrained);
plotRelativeResiduals(nexttile(tl, 4), S.relativeCurves, S.relativeNames);
plotModelComparison(nexttile(tl, 5), F);
plotTextSummary(nexttile(tl, 6), S);

end

%% ------------------------------------------------------------------------
function plotChoiceCurve(ax, curve)

hold(ax, 'on');
plot(ax, curve.meanNoise, curve.pObserved, 'ko', ...
  'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
plot(ax, curve.meanNoise, curve.pMainBinned, '-', ...
  'Color', [0.4 0.4 0.4], 'LineWidth', 1.3, 'DisplayName', 'Main');
plot(ax, curve.meanNoise, curve.pOpponentBinned, '--', ...
  'Color', [0.15 0.60 0.35], 'LineWidth', 1.4, 'DisplayName', 'Ratio 2.5');
plot(ax, curve.meanNoise, curve.pSignedBinned, '-', ...
  'Color', [0.75 0.25 0.65], 'LineWidth', 1.3, 'DisplayName', 'Free signed');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'No-change stream predictor');
ylabel(ax, 'P(correct)');
title(ax, 'All Three Streams');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotResidualCurve(ax, curve)

hold(ax, 'on');
plot(ax, curve.meanNoise, curve.residualMain, 'o-', ...
  'Color', [0.35 0.35 0.35], 'MarkerFaceColor', [0.35 0.35 0.35], ...
  'LineWidth', 1.1, 'DisplayName', 'Main residual');
plot(ax, curve.meanNoise, curve.residualOpponent, 'o--', ...
  'Color', [0.15 0.60 0.35], 'MarkerFaceColor', [0.15 0.60 0.35], ...
  'LineWidth', 1.1, 'DisplayName', 'Ratio 2.5 residual');
plot(ax, curve.meanNoise, curve.residualSigned, 'o-', ...
  'Color', [0.75 0.25 0.65], 'MarkerFaceColor', [0.75 0.25 0.65], ...
  'LineWidth', 1.1, 'DisplayName', 'Free signed residual');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'No-change stream predictor');
ylabel(ax, 'Observed - predicted');
title(ax, 'All-Stream Residuals');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotSignedSlopes(ax, slopes, opponentFit)

x = 1:height(slopes);
freeGain = slopes.gain;
constrainedGain = [opponentFit.gNegative; opponentFit.gPositive];
values = [freeGain, constrainedGain];
bar(ax, x, values, 'grouped');
hold(ax, 'on');
low = freeGain - slopes.CI95Low;
high = slopes.CI95High - freeGain;
errorbar(ax, x - 0.15, freeGain, low, high, 'k.', 'LineWidth', 1.3);
yline(ax, 0, 'k:');
set(ax, 'XTick', x, 'XTickLabel', slopes.signName);
xlim(ax, [0.25 height(slopes) + 0.75]);
ylabel(ax, 'Noise gain');
title(ax, 'No-Change Signed Gains');
legend(ax, {'Free', 'Ratio 2.5'}, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotRelativeResiduals(ax, curves, names)

colors = lines(3);
hold(ax, 'on');
for iDir = 1:3
  plot(ax, curves{iDir}.meanNoise, curves{iDir}.residualMain, 'o-', ...
    'Color', colors(iDir, :), 'MarkerFaceColor', colors(iDir, :), ...
    'LineWidth', 1.0, 'DisplayName', names{iDir});
end
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'No-change stream predictor');
ylabel(ax, 'Observed - main prediction');
title(ax, 'Relative-Direction Residuals');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotModelComparison(ax, F)

fits = {F.main, F.opponentConstrained, F.signedAsymmetry};
labels = {'Linear', 'Ratio 2.5', 'Free signed'};
deltaAIC = cellfun(@(x) x.deltaAIC, fits);
bar(ax, 1:3, deltaAIC, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'k');
set(ax, 'XTick', 1:3, 'XTickLabel', labels);
ylabel(ax, '\DeltaAIC from best');
title(ax, 'Opponent Model Comparison');
grid(ax, 'on');
box(ax, 'off');
text(ax, 2, deltaAIC(2), sprintf('vs free p=%.3g', ...
  F.opponentConstrained.pValueLRTVsFree), ...
  'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
  'FontSize', 8);

end

%% ------------------------------------------------------------------------
function plotTextSummary(ax, S)

F = S.signedFits;
slopes = F.signedSlopes;
opponent = F.opponentConstrained;
freeRatio = F.freePreferredToNullRatio;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s; side: %s', S.noisePredictorType, S.side)
  'Physical coherence signs retained'
  ''
  sprintf('Main g_{All}: %.3f [%.3f, %.3f]', F.main.gain, ...
    F.main.CI95(1, 1), F.main.CI95(1, 2))
  sprintf('Asymmetry h_{All}: %.3f [%.3f, %.3f]', ...
    F.signedAsymmetry.gain(2), F.signedAsymmetry.CI95(2, 1), ...
    F.signedAsymmetry.CI95(2, 2))
  sprintf('g_- %.3f; g_+ %.3f', slopes.gain(1), slopes.gain(2))
  sprintf('\\DeltaNLL %.3f; LRT p %.3g', ...
    F.signedAsymmetry.deltaNLL, F.signedAsymmetry.pValueLRT)
  ''
  'Opponent ratio 2.5:'
  sprintf('  k_P %.3f; k_N %.3f', opponent.kPreferred, opponent.kNull)
  sprintf('  NLL %.3f; vs free p %.3g', ...
    opponent.negLogLikelihood, opponent.pValueLRTVsFree)
  sprintf('Free preferred:null ratio %.3f', freeRatio.estimate)
  sprintf('  95%% CI [%.3f, %.3f]', freeRatio.CI95(1), freeRatio.CI95(2))
  ''
  'Max prediction:'
  '  g_+ more negative than g_-'
  '  equivalently, h_{All} < 0'
  };
text(ax, -0.05, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7.5);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.signedFits;
slopes = F.signedSlopes;
opponent = F.opponentConstrained;
freeRatio = F.freePreferredToNullRatio;
fprintf('IDQ no-change-side signed-noise analysis (%s predictor)\n', S.noisePredictorType);
fprintf('Main: gAll %.6f, NLL %.6f\n', ...
  F.main.gain, F.main.negLogLikelihood);
fprintf(['Signed asymmetry: hAll %.6f, gNegative %.6f, gPositive %.6f, ' ...
  'deltaNLL %.6f, LRT p %.6g\n'], ...
  F.signedAsymmetry.gain(2), slopes.gain(1), slopes.gain(2), ...
  F.signedAsymmetry.deltaNLL, F.signedAsymmetry.pValueLRT);
fprintf(['Opponent ratio 2.5: kPreferred %.6f, kNull %.6f, ' ...
  'gNegative %.6f, gPositive %.6f, NLL %.6f, ' ...
  'deltaNLL vs free %.6f, LRT p %.6g, deltaAIC %.6f\n'], ...
  opponent.kPreferred, opponent.kNull, opponent.gNegative, ...
  opponent.gPositive, opponent.negLogLikelihood, ...
  opponent.deltaNLLToFree, opponent.pValueLRTVsFree, opponent.deltaAIC);
fprintf('Free preferred:null ratio %.6f, SE %.6f, 95%% CI [%.6f, %.6f]\n', ...
  freeRatio.estimate, freeRatio.SE, freeRatio.CI95(1), freeRatio.CI95(2));

end
