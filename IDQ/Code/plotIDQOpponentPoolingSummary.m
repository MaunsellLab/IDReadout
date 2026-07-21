function plotIDQOpponentPoolingSummary()
% plotIDQOpponentPoolingSummary
%
% Compare opponent-rectified summation and max pooling in IDQ. The primary
% test isolates no-change pooling while fitting the empirical change-side
% drift and pooled non-drift gains. A secondary 2-by-2 comparison applies
% sum/max pooling independently to the two sides.
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_OpponentPoolingSummary.mat
%   Plots/AcrossSessionSummaries/IDQOpponentPoolingSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQOpponentPoolingSummary:NoSessionAnalysisFiles', ...
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
preferredToNullRatio = 2.42;

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

poolingFits = fitIDQOpponentPooling(trialTable, psychDir.sessionFits, ...
  alignedWeibull, targetPerformance, preferredToNullRatio);
disagreementCurve = makeDisagreementCurve(poolingFits.trialData, 9);

opponentPoolingSummary = struct();
opponentPoolingSummary.createdBy = mfilename;
opponentPoolingSummary.createdAt = datetime('now');
opponentPoolingSummary.noisePredictorType = 'rectStep';
opponentPoolingSummary.nSessions = numel(sessionAnalyses);
opponentPoolingSummary.nTrials = height(trialTable);
opponentPoolingSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
opponentPoolingSummary.sessionAnalysisFolder = sessionAnalysisFolder;
opponentPoolingSummary.initialBetaWeibull = initialBetaWeibull;
opponentPoolingSummary.initialLapse = initialLapse;
opponentPoolingSummary.targetPerformance = targetPerformance;
opponentPoolingSummary.lapseBounds = lapseBounds;
opponentPoolingSummary.preferredToNullRatio = preferredToNullRatio;
opponentPoolingSummary.sessionFits = psychDir.sessionFits;
opponentPoolingSummary.alignedWeibull = alignedWeibull;
opponentPoolingSummary.trialTable = trialTable;
opponentPoolingSummary.poolingFits = poolingFits;
opponentPoolingSummary.disagreementCurve = disagreementCurve;

summaryFile = fullfile(summaryFolder, 'IDQ_OpponentPoolingSummary.mat');
save(summaryFile, 'opponentPoolingSummary', '-v7.3');

fig = plotPoolingFigure(opponentPoolingSummary);
pdfFile = fullfile(plotFolder, 'IDQOpponentPoolingSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
printFitSummary(opponentPoolingSummary);

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
    error('plotIDQOpponentPoolingSummary:OldSessionAnalysis', ...
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
function curve = makeDisagreementCurve(T, nBins)

zSum = standardize(T.noChangeSum);
zMax = standardize(T.noChangeMax);
disagreement = zMax - zSum;
edges = quantile(disagreement, linspace(0, 1, nBins + 1));
edges(1) = -Inf;
edges(end) = Inf;
bin = discretize(disagreement, edges);

meanDisagreement = nan(nBins, 1);
residualSum = nan(nBins, 1);
residualMax = nan(nBins, 1);
pObserved = nan(nBins, 1);
pSum = nan(nBins, 1);
pMax = nan(nBins, 1);
nTrials = zeros(nBins, 1);
for iBin = 1:nBins
  idx = bin == iBin;
  nTrials(iBin) = sum(idx);
  meanDisagreement(iBin) = mean(disagreement(idx), 'omitnan');
  residualSum(iBin) = mean(T.correct(idx) - T.pControlledSum(idx), 'omitnan');
  residualMax(iBin) = mean(T.correct(idx) - T.pControlledMax(idx), 'omitnan');
  pObserved(iBin) = mean(T.correct(idx), 'omitnan');
  pSum(iBin) = mean(T.pControlledSum(idx), 'omitnan');
  pMax(iBin) = mean(T.pControlledMax(idx), 'omitnan');
end
curve = table(meanDisagreement, residualSum, residualMax, ...
  pObserved, pSum, pMax, nTrials);

end

%% ------------------------------------------------------------------------
function z = standardize(x)

sigma = std(x, 0, 'omitnan');
if ~isfinite(sigma) || sigma <= 0
  z = zeros(size(x));
else
  z = (x - mean(x, 'omitnan')) / sigma;
end

end

%% ------------------------------------------------------------------------
function fig = plotPoolingFigure(S)

F = S.poolingFits;
fig = figure(550);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf(['IDQ Opponent-Rectified Sum versus Max Pooling ' ...
  '(%d sessions; preferred:null %.2f)'], S.nSessions, ...
  S.preferredToNullRatio), 'Interpreter', 'none', 'FontWeight', 'bold');

plotControlledComparison(nexttile(tl, 1), F);
plotControlledGains(nexttile(tl, 2), F);
plotDisagreement(nexttile(tl, 3), S.disagreementCurve);
plotTwoByTwo(nexttile(tl, 4), F);
plotWinnerFractions(nexttile(tl, 5), F.winnerStats);
plotTextSummary(nexttile(tl, 6), S);

end

%% ------------------------------------------------------------------------
function plotControlledComparison(ax, F)

fits = {F.empiricalLinear, F.empiricalFreeSigned, ...
  F.controlledSum, F.controlledMax};
labels = {'Linear', 'Free signed', 'Opponent sum', 'Opponent max'};
deltaAIC = cellfun(@(x) x.deltaAIC, fits);
bar(ax, 1:numel(fits), deltaAIC, ...
  'FaceColor', [0.50 0.62 0.74], 'EdgeColor', 'k');
set(ax, 'XTick', 1:numel(fits), 'XTickLabel', labels);
xtickangle(ax, 25);
ylabel(ax, '\DeltaAIC from best');
title(ax, 'Controlled No-Change Test');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotControlledGains(ax, F)

values = [F.controlledSum.gain, F.controlledMax.gain];
b = bar(ax, 1:3, values, 'grouped');
b(1).FaceColor = [0.35 0.60 0.82];
b(2).FaceColor = [0.82 0.45 0.35];
hold(ax, 'on');
for iModel = 1:2
  fit = {F.controlledSum, F.controlledMax};
  thisFit = fit{iModel};
  xOffset = (iModel - 1.5) * 0.28;
  errorbar(ax, (1:3) + xOffset, thisFit.gain, thisFit.SE, ...
    'k.', 'LineWidth', 1.1);
end
yline(ax, 0, 'k:');
set(ax, 'XTick', 1:3, 'XTickLabel', {'Drift', 'Change non-drift', 'NC pool'});
ylabel(ax, 'Fitted gain');
title(ax, 'Equal-Parameter Fits');
legend(ax, {'Sum', 'Max'}, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotDisagreement(ax, curve)

hold(ax, 'on');
plot(ax, curve.meanDisagreement, curve.residualSum, 'o-', ...
  'Color', [0.35 0.60 0.82], 'MarkerFaceColor', [0.35 0.60 0.82], ...
  'LineWidth', 1.1, 'DisplayName', 'Sum residual');
plot(ax, curve.meanDisagreement, curve.residualMax, 'o-', ...
  'Color', [0.82 0.45 0.35], 'MarkerFaceColor', [0.82 0.45 0.35], ...
  'LineWidth', 1.1, 'DisplayName', 'Max residual');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'z(max pool) - z(sum pool)');
ylabel(ax, 'Observed - predicted');
title(ax, 'No-Change Pooling Disagreement');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotTwoByTwo(ax, F)

order = F.pooling2x2Order;
fits = cell(size(order));
for iModel = 1:numel(order)
  fits{iModel} = F.pooling2x2.(order{iModel});
end
labels = {'C sum / N sum', 'C sum / N max', ...
  'C max / N sum', 'C max / N max'};
deltaAIC = cellfun(@(x) x.deltaAIC, fits);
bar(ax, 1:numel(fits), deltaAIC, ...
  'FaceColor', [0.62 0.55 0.72], 'EdgeColor', 'k');
set(ax, 'XTick', 1:numel(fits), 'XTickLabel', labels);
xtickangle(ax, 25);
ylabel(ax, '\DeltaAIC from best');
title(ax, 'Two-Side 2 x 2 Pooling Test');
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotWinnerFractions(ax, W)

values = [W.changeDriftMatchedFraction, W.changeNonDriftFraction; ...
  W.noChangePreferredFraction, W.noChangeOpponentFraction];
bar(ax, values, 'stacked');
set(ax, 'XTick', 1:2, 'XTickLabel', {'Change', 'No-change'});
ylim(ax, [0 1]);
ylabel(ax, 'Fraction of trials');
title(ax, 'Latent Max Candidate');
legend(ax, {'Drift / preferred', 'Non-drift / opponent'}, ...
  'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotTextSummary(ax, S)

F = S.poolingFits;
C = F.controlledContrast;
P = F.pooling2x2;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s; preferred:null %.2f', ...
    S.noisePredictorType, S.preferredToNullRatio)
  'Pools centered within session'
  ''
  'Controlled no-change comparison:'
  sprintf('  Sum NLL %.3f; AIC %.3f', ...
    F.controlledSum.negLogLikelihood, F.controlledSum.AIC)
  sprintf('  Max NLL %.3f; AIC %.3f', ...
    F.controlledMax.negLogLikelihood, F.controlledMax.AIC)
  sprintf('  DeltaNLL(max over sum) %.3f', C.deltaNLLMaxOverSum)
  sprintf('  Sessions favoring max %d/%d', ...
    C.nSessionsFavoringMax, numel(C.sessionIndex))
  sprintf('  Corr(sum,max) %.3f', F.poolCorrelation(3, 4))
  ''
  'Two-side comparison:'
  sprintf('  Csum/Nsum NLL %.3f', P.changeSum_noChangeSum.negLogLikelihood)
  sprintf('  Csum/Nmax NLL %.3f', P.changeSum_noChangeMax.negLogLikelihood)
  sprintf('  Cmax/Nsum NLL %.3f', P.changeMax_noChangeSum.negLogLikelihood)
  sprintf('  Cmax/Nmax NLL %.3f', P.changeMax_noChangeMax.negLogLikelihood)
  ''
  sprintf('Change max non-drift winner %.4f', ...
    F.winnerStats.changeNonDriftFraction)
  'Winner labels are model-implied, not observed'
  };
text(ax, -0.05, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7.2);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.poolingFits;
C = F.controlledContrast;
fprintf(['IDQ opponent pooling analysis (%s predictor; ' ...
  'preferred:null %.2f)\n'], S.noisePredictorType, S.preferredToNullRatio);
fprintf(['Controlled empirical change + no-change sum: ' ...
  'gD %.6f, gChangeN %.6f, gPool %.6f, NLL %.6f, AIC %.6f\n'], ...
  F.controlledSum.gain(1), F.controlledSum.gain(2), ...
  F.controlledSum.gain(3), F.controlledSum.negLogLikelihood, ...
  F.controlledSum.AIC);
fprintf(['Controlled empirical change + no-change max: ' ...
  'gD %.6f, gChangeN %.6f, gPool %.6f, NLL %.6f, AIC %.6f\n'], ...
  F.controlledMax.gain(1), F.controlledMax.gain(2), ...
  F.controlledMax.gain(3), F.controlledMax.negLogLikelihood, ...
  F.controlledMax.AIC);
fprintf(['Controlled max over sum: deltaNLL %.6f; ' ...
  'sessions favoring max %d/%d; corr(sum,max) %.6f\n'], ...
  C.deltaNLLMaxOverSum, C.nSessionsFavoringMax, numel(C.sessionIndex), ...
  F.poolCorrelation(3, 4));
fprintf('Empirical linear: NLL %.6f, deltaAIC %.6f\n', ...
  F.empiricalLinear.negLogLikelihood, F.empiricalLinear.deltaAIC);
fprintf('Empirical free signed: NLL %.6f, deltaAIC %.6f\n', ...
  F.empiricalFreeSigned.negLogLikelihood, F.empiricalFreeSigned.deltaAIC);

order = F.pooling2x2Order;
for iModel = 1:numel(order)
  fit = F.pooling2x2.(order{iModel});
  fprintf(['Two-side %s: gChange %.6f, gNoChange %.6f, ' ...
    'NLL %.6f, deltaAIC %.6f\n'], order{iModel}, fit.gain(1), ...
    fit.gain(2), fit.negLogLikelihood, fit.deltaAIC);
end
fprintf(['Two-side pooling contrasts: max over sum on change %.6f; ' ...
  'max over sum on no-change %.6f\n'], ...
  F.changePoolingContrast.deltaNLLMaxOverSum, ...
  F.noChangePoolingContrast.deltaNLLMaxOverSum);
fprintf(['Latent max winners: change drift-matched %.6f, ' ...
  'change non-drift %.6f; no-change preferred %.6f, opponent %.6f\n'], ...
  F.winnerStats.changeDriftMatchedFraction, ...
  F.winnerStats.changeNonDriftFraction, ...
  F.winnerStats.noChangePreferredFraction, ...
  F.winnerStats.noChangeOpponentFraction);

end
