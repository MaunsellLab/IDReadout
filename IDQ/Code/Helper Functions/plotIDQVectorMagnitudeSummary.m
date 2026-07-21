function plotIDQVectorMagnitudeSummary()
% plotIDQVectorMagnitudeSummary
%
% Direct test of a net-motion-vector-magnitude strategy for IDQ choices.
% Uses the six step-rectangular stream predictors but does not fit the full
% six-gain additive model.
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_VectorMagnitudeSummary.mat
%   Plots/AcrossSessionSummaries/IDQVectorMagnitudeSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQVectorMagnitudeSummary:NoSessionAnalysisFiles', ...
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
  trialTable.(sprintf('changeNoisePredDir%d', iDir)) = ...
    trialTable.(sprintf('rectStepChangeNoisePredDir%d', iDir));
  trialTable.(sprintf('noChangeNoisePredDir%d', iDir)) = ...
    trialTable.(sprintf('rectStepNoChangeNoisePredDir%d', iDir));
end

vectorFits = fitIDQVectorMagnitude(trialTable, ...
  psychDir.sessionFits, alignedWeibull, targetPerformance);
vectorCurve = makeBinnedCurve(vectorFits.trialData, 'xSharedVector', 10);
changeCurve = makeBinnedCurve(vectorFits.trialData, 'xChangeVector', 10);
noChangeCurve = makeBinnedCurve(vectorFits.trialData, 'rNoChange', 10);

vectorMagnitudeSummary = struct();
vectorMagnitudeSummary.createdBy = mfilename;
vectorMagnitudeSummary.createdAt = datetime('now');
vectorMagnitudeSummary.noisePredictorType = 'rectStep';
vectorMagnitudeSummary.nSessions = numel(sessionAnalyses);
vectorMagnitudeSummary.nTrials = height(trialTable);
vectorMagnitudeSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
vectorMagnitudeSummary.sessionAnalysisFolder = sessionAnalysisFolder;
vectorMagnitudeSummary.initialBetaWeibull = initialBetaWeibull;
vectorMagnitudeSummary.initialLapse = initialLapse;
vectorMagnitudeSummary.targetPerformance = targetPerformance;
vectorMagnitudeSummary.lapseBounds = lapseBounds;
vectorMagnitudeSummary.sessionFits = psychDir.sessionFits;
vectorMagnitudeSummary.alignedWeibull = alignedWeibull;
vectorMagnitudeSummary.trialTable = trialTable;
vectorMagnitudeSummary.vectorFits = vectorFits;
vectorMagnitudeSummary.vectorCurve = vectorCurve;
vectorMagnitudeSummary.changeCurve = changeCurve;
vectorMagnitudeSummary.noChangeCurve = noChangeCurve;

summaryFile = fullfile(summaryFolder, 'IDQ_VectorMagnitudeSummary.mat');
save(summaryFile, 'vectorMagnitudeSummary', '-v7.3');

fig = plotVectorFigure(vectorMagnitudeSummary);
pdfFile = fullfile(plotFolder, 'IDQVectorMagnitudeSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
printFitSummary(vectorMagnitudeSummary);

end

%% ------------------------------------------------------------------------
function validateInputs(sessionAnalyses)

requiredVars = strings(1, 6);
iVar = 0;
for sideToken = {'Change', 'NoChange'}
  for iDir = 1:3
    iVar = iVar + 1;
    requiredVars(iVar) = sprintf('rectStep%sNoisePredDir%d', ...
      sideToken{1}, iDir);
  end
end
for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  missing = requiredVars(~ismember(requiredVars, ...
    string(SA.trialTable.Properties.VariableNames)));
  if ~isempty(missing)
    error('plotIDQVectorMagnitudeSummary:OldSessionAnalysis', ...
      ['Session %s lacks the six side-specific rectangular predictors. ' ...
       'Rerun makeIDQSessionAnalyses before this analysis.'], SA.fileName);
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
function curve = makeBinnedCurve(T, predictorName, nBins)

x = T.(predictorName);
edges = quantile(x, linspace(0, 1, nBins + 1));
edges(1) = -Inf;
edges(end) = Inf;
bin = discretize(x, edges);

meanPredictor = nan(nBins, 1);
pObserved = nan(nBins, 1);
pStepOnly = nan(nBins, 1);
pVectorFixed = nan(nBins, 1);
pVectorShared = nan(nBins, 1);
pVectorSide = nan(nBins, 1);
pAdditiveReduced = nan(nBins, 1);
pAdditiveSigned = nan(nBins, 1);
residualVectorSide = nan(nBins, 1);
residualAdditiveReduced = nan(nBins, 1);
residualAdditiveSigned = nan(nBins, 1);
nTrials = zeros(nBins, 1);
for iBin = 1:nBins
  idx = bin == iBin;
  nTrials(iBin) = sum(idx);
  meanPredictor(iBin) = mean(x(idx), 'omitnan');
  pObserved(iBin) = mean(T.correct(idx), 'omitnan');
  pStepOnly(iBin) = mean(T.pStepOnly(idx), 'omitnan');
  pVectorFixed(iBin) = mean(T.pVectorFixed(idx), 'omitnan');
  pVectorShared(iBin) = mean(T.pVectorShared(idx), 'omitnan');
  pVectorSide(iBin) = mean(T.pVectorSide(idx), 'omitnan');
  pAdditiveReduced(iBin) = mean(T.pAdditiveReduced(idx), 'omitnan');
  pAdditiveSigned(iBin) = mean(T.pAdditiveSigned(idx), 'omitnan');
  residualVectorSide(iBin) = mean(T.correct(idx) - T.pVectorSide(idx), 'omitnan');
  residualAdditiveReduced(iBin) = mean(T.correct(idx) - T.pAdditiveReduced(idx), 'omitnan');
  residualAdditiveSigned(iBin) = mean(T.correct(idx) - T.pAdditiveSigned(idx), 'omitnan');
end
curve = table(meanPredictor, pObserved, pStepOnly, pVectorFixed, ...
  pVectorShared, pVectorSide, pAdditiveReduced, pAdditiveSigned, ...
  residualVectorSide, residualAdditiveReduced, residualAdditiveSigned, nTrials);
curve.Properties.Description = predictorName;

end

%% ------------------------------------------------------------------------
function fig = plotVectorFigure(S)

F = S.vectorFits;
fig = figure(540);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Net Motion-Vector Magnitude (%d sessions)', S.nSessions), ...
  'Interpreter', 'none', 'FontWeight', 'bold');

plotNLLImprovement(nexttile(tl, 1), F);
plotDeltaAIC(nexttile(tl, 2), F);
plotVectorGains(nexttile(tl, 3), F.vectorSideGains);
plotChoiceCurve(nexttile(tl, 4), S.vectorCurve);
plotResidualCurve(nexttile(tl, 5), S.changeCurve, ...
  'Change-Vector Perturbation', 'Change vector perturbation');
plotResidualCurve(nexttile(tl, 6), S.noChangeCurve, ...
  'No-Change Vector Magnitude', 'No-change vector magnitude');
plotImpliedChangeGains(nexttile(tl, 7), F);
plotMagnitudeDistributions(nexttile(tl, 8), F.trialData);
plotTextSummary(nexttile(tl, 9), S);

end

%% ------------------------------------------------------------------------
function [fits, labels] = orderedFits(F)

fits = {F.stepOnly, F.vectorFixed, F.vectorSharedGain, ...
  F.vectorSideGains, F.additiveReduced, F.additiveSigned};
labels = {'Step', 'Vector fixed', 'Vector shared', 'Vector sides', ...
  'Additive', 'Add.+signed'};

end

%% ------------------------------------------------------------------------
function plotNLLImprovement(ax, F)

[fits, labels] = orderedFits(F);
stepNLL = F.stepOnly.negLogLikelihood;
improvement = cellfun(@(x) stepNLL - x.negLogLikelihood, fits);
bar(ax, 1:numel(fits), improvement, 'FaceColor', [0.45 0.60 0.75], 'EdgeColor', 'k');
set(ax, 'XTick', 1:numel(fits), 'XTickLabel', labels);
xtickangle(ax, 30);
ylabel(ax, 'NLL improvement over step only');
title(ax, 'Predictive Improvement');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotDeltaAIC(ax, F)

[fits, labels] = orderedFits(F);
deltaAIC = cellfun(@(x) x.deltaAIC, fits);
bar(ax, 1:numel(fits), deltaAIC, 'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'k');
set(ax, 'XTick', 1:numel(fits), 'XTickLabel', labels);
xtickangle(ax, 30);
ylabel(ax, '\DeltaAIC from best model');
title(ax, 'Model Comparison');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotVectorGains(ax, fit)

x = 1:2;
y = fit.gain;
low = y - fit.CI95(:, 1);
high = fit.CI95(:, 2) - y;
bar(ax, x, y, 'BarWidth', 0.55, 'FaceColor', [0.45 0.70 0.55], 'EdgeColor', 'k');
hold(ax, 'on');
errorbar(ax, x, y, low, high, 'k.', 'LineWidth', 1.3);
yline(ax, 1, 'k:', 'Physical-vector value');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
set(ax, 'XTick', x, 'XTickLabel', fit.predictorNames);
xlim(ax, [0.25 2.75]);
ylabel(ax, 'Vector gain');
title(ax, 'Separate Side Gains');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotChoiceCurve(ax, curve)

hold(ax, 'on');
plot(ax, curve.meanPredictor, curve.pObserved, 'ko', ...
  'MarkerFaceColor', 'k', 'DisplayName', 'Observed');
plot(ax, curve.meanPredictor, curve.pVectorFixed, '--', ...
  'Color', [0.25 0.55 0.35], 'LineWidth', 1.1, 'DisplayName', 'Vector fixed');
plot(ax, curve.meanPredictor, curve.pVectorSide, '-', ...
  'Color', [0.15 0.65 0.35], 'LineWidth', 1.4, 'DisplayName', 'Vector sides');
plot(ax, curve.meanPredictor, curve.pAdditiveSigned, '-', ...
  'Color', [0.55 0.30 0.75], 'LineWidth', 1.4, 'DisplayName', 'Add.+signed');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'Vector decision perturbation');
ylabel(ax, 'P(correct)');
title(ax, 'Binned Vector Decision Variable');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotResidualCurve(ax, curve, plotTitle, xLabelText)

hold(ax, 'on');
plot(ax, curve.meanPredictor, curve.residualVectorSide, 'o-', ...
  'Color', [0.15 0.65 0.35], 'MarkerFaceColor', [0.15 0.65 0.35], ...
  'LineWidth', 1.1, 'DisplayName', 'Vector sides');
plot(ax, curve.meanPredictor, curve.residualAdditiveSigned, 'o-', ...
  'Color', [0.55 0.30 0.75], 'MarkerFaceColor', [0.55 0.30 0.75], ...
  'LineWidth', 1.1, 'DisplayName', 'Add.+signed');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, xLabelText);
ylabel(ax, 'Observed - predicted');
title(ax, plotTitle);
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotImpliedChangeGains(ax, F)

gChange = F.vectorSideGains.gain(1);
vectorValues = [gChange, -0.5 * gChange];
additiveValues = F.additiveSigned.gain(1:2)';
values = [vectorValues; additiveValues];
bar(ax, values', 'grouped');
yline(ax, 0, 'k:');
set(ax, 'XTick', 1:2, 'XTickLabel', {'Drift', 'Each non-drift'});
ylabel(ax, 'Equivalent change-side gain');
title(ax, 'Geometry Versus Additive Fit');
legend(ax, {'Vector implied', 'Additive fitted'}, ...
  'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotMagnitudeDistributions(ax, T)

hold(ax, 'on');
histogram(ax, T.xChangeVector, 30, 'Normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'LineWidth', 1.3, 'DisplayName', 'R_C - step');
histogram(ax, T.rNoChange, 30, 'Normalization', 'probability', ...
  'DisplayStyle', 'stairs', 'LineWidth', 1.3, 'DisplayName', 'R_N');
xlabel(ax, 'Coherence-equivalent magnitude');
ylabel(ax, 'Probability');
title(ax, 'Vector Components');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotTextSummary(ax, S)

F = S.vectorFits;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s', S.noisePredictorType)
  ''
  sprintf('Step-only NLL: %.3f', F.stepOnly.negLogLikelihood)
  sprintf('Fixed-vector NLL: %.3f', F.vectorFixed.negLogLikelihood)
  sprintf('Shared-vector: g %.3f; NLL %.3f', ...
    F.vectorSharedGain.gain, F.vectorSharedGain.negLogLikelihood)
  sprintf('Side-vector: g_C %.3f, g_N %.3f', ...
    F.vectorSideGains.gain(1), F.vectorSideGains.gain(2))
  sprintf('  NLL %.3f; vs shared p %.3g', ...
    F.vectorSideGains.negLogLikelihood, F.vectorSideGains.pValueLRT)
  ''
  sprintf('Add.+signed: g_D %.3f, g_CN %.3f', ...
    F.additiveSigned.gain(1), F.additiveSigned.gain(2))
  sprintf('  g_N %.3f, h_N %.3f; NLL %.3f', ...
    F.additiveSigned.gain(3), F.additiveSigned.gain(4), ...
    F.additiveSigned.negLogLikelihood)
  ''
  'R^2=a^2+b^2+c^2-ab-ac-bc'
  };
text(ax, -0.10, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.vectorFits;
fprintf('IDQ vector-magnitude analysis (%s predictor)\n', S.noisePredictorType);
fprintf('Step only: NLL %.6f, AIC %.6f, deltaAIC %.6f\n', ...
  F.stepOnly.negLogLikelihood, F.stepOnly.AIC, F.stepOnly.deltaAIC);
fprintf('Vector fixed: NLL %.6f, AIC %.6f, deltaAIC %.6f\n', ...
  F.vectorFixed.negLogLikelihood, F.vectorFixed.AIC, F.vectorFixed.deltaAIC);
fprintf(['Vector shared gain: g %.6f, NLL %.6f, deltaNLL vs fixed %.6f, ' ...
  'LRT p %.6g, deltaAIC %.6f\n'], F.vectorSharedGain.gain, ...
  F.vectorSharedGain.negLogLikelihood, F.vectorSharedGain.deltaNLL, ...
  F.vectorSharedGain.pValueLRT, F.vectorSharedGain.deltaAIC);
fprintf(['Vector side gains: gChange %.6f, gNoChange %.6f, NLL %.6f, ' ...
  'deltaNLL vs shared %.6f, LRT p %.6g, deltaAIC %.6f\n'], ...
  F.vectorSideGains.gain(1), F.vectorSideGains.gain(2), ...
  F.vectorSideGains.negLogLikelihood, F.vectorSideGains.deltaNLL, ...
  F.vectorSideGains.pValueLRT, F.vectorSideGains.deltaAIC);
fprintf(['Additive reduced: gD %.6f, gChangeN %.6f, gNoChange %.6f, ' ...
  'NLL %.6f, deltaAIC %.6f\n'], F.additiveReduced.gain(1), ...
  F.additiveReduced.gain(2), F.additiveReduced.gain(3), ...
  F.additiveReduced.negLogLikelihood, F.additiveReduced.deltaAIC);
fprintf(['Additive + signed no-change: gD %.6f, gChangeN %.6f, ' ...
  'gNoChange %.6f, hNoChange %.6f, NLL %.6f, deltaNLL vs additive %.6f, ' ...
  'LRT p %.6g, deltaAIC %.6f\n'], F.additiveSigned.gain(1), ...
  F.additiveSigned.gain(2), F.additiveSigned.gain(3), ...
  F.additiveSigned.gain(4), F.additiveSigned.negLogLikelihood, ...
  F.additiveSigned.deltaNLL, F.additiveSigned.pValueLRT, ...
  F.additiveSigned.deltaAIC);
fprintf('Vector side model implied change gains: drift %.6f, each non-drift %.6f\n', ...
  F.vectorSideGains.gain(1), -0.5 * F.vectorSideGains.gain(1));

end
