function plotIDQOpponentPNormSummary()
% plotIDQOpponentPNormSummary
%
% Fit a one-parameter p-norm family spanning opponent-rectified summation
% (p=1) to hard max (p=Inf) on the no-change side. The primary mechanistic
% model holds change-side pooling at hard max. A controlled model instead
% retains empirical change drift and pooled non-drift terms.
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_OpponentPNormSummary.mat
%   Plots/AcrossSessionSummaries/IDQOpponentPNormSummary.pdf

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
  error('plotIDQOpponentPNormSummary:NoSessionAnalysisFiles', ...
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

pNormFits = fitIDQOpponentPNorm(trialTable, psychDir.sessionFits, ...
  alignedWeibull, targetPerformance, preferredToNullRatio);

opponentPNormSummary = struct();
opponentPNormSummary.createdBy = mfilename;
opponentPNormSummary.createdAt = datetime('now');
opponentPNormSummary.noisePredictorType = 'rectStep';
opponentPNormSummary.nSessions = numel(sessionAnalyses);
opponentPNormSummary.nTrials = height(trialTable);
opponentPNormSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
opponentPNormSummary.sessionAnalysisFolder = sessionAnalysisFolder;
opponentPNormSummary.initialBetaWeibull = initialBetaWeibull;
opponentPNormSummary.initialLapse = initialLapse;
opponentPNormSummary.targetPerformance = targetPerformance;
opponentPNormSummary.lapseBounds = lapseBounds;
opponentPNormSummary.preferredToNullRatio = preferredToNullRatio;
opponentPNormSummary.sessionFits = psychDir.sessionFits;
opponentPNormSummary.alignedWeibull = alignedWeibull;
opponentPNormSummary.trialTable = trialTable;
opponentPNormSummary.pNormFits = pNormFits;

summaryFile = fullfile(summaryFolder, 'IDQ_OpponentPNormSummary.mat');
save(summaryFile, 'opponentPNormSummary', '-v7.3');

fig = plotPNormFigure(opponentPNormSummary);
pdfFile = fullfile(plotFolder, 'IDQOpponentPNormSummary.pdf');
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
printFitSummary(opponentPNormSummary);

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
    error('plotIDQOpponentPNormSummary:OldSessionAnalysis', ...
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
function fig = plotPNormFigure(S)

F = S.pNormFits;
fig = figure(560);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf(['IDQ Opponent-Rectified p-Norm Pooling ' ...
  '(%d sessions; preferred:null %.2f)'], S.nSessions, ...
  S.preferredToNullRatio), 'Interpreter', 'none', 'FontWeight', 'bold');

plotLikelihoodProfile(nexttile(tl, 1), F.mechanistic, ...
  'Change Max + No-Change p-Norm');
plotLikelihoodProfile(nexttile(tl, 2), F.controlled, ...
  'Empirical Change + No-Change p-Norm');
plotAICComparison(nexttile(tl, 3), F);
plotGainProfile(nexttile(tl, 4), F);
plotSessionContrasts(nexttile(tl, 5), F);
plotTextSummary(nexttile(tl, 6), S);

end


%% ------------------------------------------------------------------------
function plotLikelihoodProfile(ax, context, titleText)

P = context.profile;
C = context.continuous;
hold(ax, 'on');
semilogx(ax, P.p, P.deltaNLL, 'o-', ...
  'Color', [0.35 0.55 0.78], 'MarkerFaceColor', [0.35 0.55 0.78], ...
  'LineWidth', 1.2, 'DisplayName', 'Fixed p');
semilogx(ax, C.p, C.negLogLikelihood - ...
  min(C.negLogLikelihood, context.hardMax.negLogLikelihood), 'p', ...
  'Color', [0.80 0.35 0.25], 'MarkerFaceColor', [0.80 0.35 0.25], ...
  'MarkerSize', 10, 'DisplayName', 'Fitted p');
xHard = max(P.p) * 1.6;
semilogx(ax, xHard, P.hardMaxDeltaNLL, 'ks', ...
  'MarkerFaceColor', 'k', 'DisplayName', 'Hard max');
yline(ax, 0.5 * 3.84145882069412, 'k:', ...
  'DisplayName', '95% profile cutoff');
set(ax, 'XTick', [1 2 4 8 16 32 64 128 xHard], ...
  'XTickLabel', {'1','2','4','8','16','32','64','128','Inf'});
xlim(ax, [1 xHard * 1.08]);
xlabel(ax, 'Pooling exponent p');
ylabel(ax, 'Profile \DeltaNLL');
title(ax, titleText);
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotAICComparison(ax, F)

values = [F.mechanistic.continuous.sumDeltaAIC, ...
          F.controlled.continuous.sumDeltaAIC; ...
          F.mechanistic.continuous.deltaAIC, ...
          F.controlled.continuous.deltaAIC; ...
          F.mechanistic.continuous.hardMaxDeltaAIC, ...
          F.controlled.continuous.hardMaxDeltaAIC];
bar(ax, values, 'grouped');
set(ax, 'XTick', 1:3, 'XTickLabel', {'Sum', 'Fitted p', 'Hard max'});
ylabel(ax, '\DeltaAIC within context');
title(ax, 'Complexity-Adjusted Comparison');
legend(ax, {'Mechanistic', 'Controlled'}, ...
  'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotGainProfile(ax, F)

hold(ax, 'on');
semilogx(ax, F.mechanistic.profile.p, F.mechanistic.profile.gain(:, end), ...
  'o-', 'Color', [0.30 0.55 0.78], 'MarkerFaceColor', [0.30 0.55 0.78], ...
  'LineWidth', 1.1, 'DisplayName', 'Mechanistic');
semilogx(ax, F.controlled.profile.p, F.controlled.profile.gain(:, end), ...
  'o-', 'Color', [0.78 0.42 0.30], 'MarkerFaceColor', [0.78 0.42 0.30], ...
  'LineWidth', 1.1, 'DisplayName', 'Controlled');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlabel(ax, 'Pooling exponent p');
ylabel(ax, 'No-change pool gain');
title(ax, 'Gain Rescaling across p');
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotSessionContrasts(ax, F)

M = F.mechanistic.endpointContrast;
C = F.controlled.endpointContrast;
values = [M.sessionDeltaNLLMaxOverSum, C.sessionDeltaNLLMaxOverSum];
bar(ax, M.sessionIndex, values, 'grouped');
yline(ax, 0, 'k:');
xlabel(ax, 'Session index');
ylabel(ax, '\DeltaNLL favoring hard max');
title(ax, 'Hard Max versus Sum by Session');
legend(ax, {'Mechanistic', 'Controlled'}, ...
  'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
box(ax, 'off');

end

%% ------------------------------------------------------------------------
function plotTextSummary(ax, S)

F = S.pNormFits;
M = F.mechanistic;
C = F.controlled;
axis(ax, 'off');
txt = {
  sprintf('Created: %s', string(S.createdAt, 'dd/MM/yyyy HH:mm'))
  sprintf('%d sessions; %d noisy trials', S.nSessions, F.nTrials)
  sprintf('Predictor: %s; preferred:null %.2f', ...
    S.noisePredictorType, S.preferredToNullRatio)
  'No-change pools centered within session'
  ''
  'Mechanistic (change hard max):'
  pEstimateText(M.continuous, M.profile)
  sprintf('  Sum NLL %.3f; hard-max NLL %.3f', ...
    M.sum.negLogLikelihood, M.hardMax.negLogLikelihood)
  sprintf('  Fitted-p NLL %.3f; AIC delta %.3f', ...
    M.continuous.negLogLikelihood, M.continuous.deltaAIC)
  sprintf('  Hard-max AIC delta %.3f', M.continuous.hardMaxDeltaAIC)
  ''
  'Controlled empirical change:'
  pEstimateText(C.continuous, C.profile)
  sprintf('  Sum NLL %.3f; hard-max NLL %.3f', ...
    C.sum.negLogLikelihood, C.hardMax.negLogLikelihood)
  sprintf('  Fitted-p NLL %.3f; AIC delta %.3f', ...
    C.continuous.negLogLikelihood, C.continuous.deltaAIC)
  sprintf('  Hard-max AIC delta %.3f', C.continuous.hardMaxDeltaAIC)
  ''
  'p=1 is sum; p=Inf is exact hard max'
  'Profile CI is grid-based; Inf means unbounded above'
  };
text(ax, -0.05, 1, txt, 'Units', 'normalized', ...
  'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontName', 'Menlo', 'FontSize', 7.2);

end

%% ------------------------------------------------------------------------
function txt = pEstimateText(fit, profile)

if fit.pAtUpperBound
  estimate = sprintf('p >= %.0f', fit.p);
else
  estimate = sprintf('p = %.3f', fit.p);
end
if isinf(profile.profileCI95High)
  ci = sprintf('[%.3g, Inf]', profile.profileCI95Low);
else
  ci = sprintf('[%.3g, %.3g]', ...
    profile.profileCI95Low, profile.profileCI95High);
end
txt = sprintf('  %s; profile 95%% CI %s', estimate, ci);

end

%% ------------------------------------------------------------------------
function printFitSummary(S)

F = S.pNormFits;
fprintf(['IDQ opponent p-norm analysis (%s predictor; ' ...
  'preferred:null %.2f)\n'], S.noisePredictorType, S.preferredToNullRatio);
printContext('Mechanistic change-max', F.mechanistic);
printContext('Controlled empirical-change', F.controlled);

end

%% ------------------------------------------------------------------------
function printContext(label, context)

C = context.continuous;
P = context.profile;
if C.pAtUpperBound
  pString = sprintf('>=%.0f', C.p);
else
  pString = sprintf('%.6f', C.p);
end
if isinf(P.profileCI95High)
  ciString = sprintf('[%.6f, Inf]', P.profileCI95Low);
else
  ciString = sprintf('[%.6f, %.6f]', ...
    P.profileCI95Low, P.profileCI95High);
end
fprintf(['%s: p %s, profile 95%% CI %s, gNoChange %.6f, ' ...
  'NLL %.6f, AIC %.6f, deltaAIC %.6f\n'], label, pString, ciString, ...
  C.gain(end), C.negLogLikelihood, C.AIC, C.deltaAIC);
fprintf(['  Sum: gNoChange %.6f, NLL %.6f, deltaAIC %.6f; ' ...
  'hard max: gNoChange %.6f, NLL %.6f, deltaAIC %.6f\n'], ...
  context.sum.gain(end), context.sum.negLogLikelihood, C.sumDeltaAIC, ...
  context.hardMax.gain(end), context.hardMax.negLogLikelihood, ...
  C.hardMaxDeltaAIC);
fprintf(['  Hard max over sum: deltaNLL %.6f; ' ...
  'sessions favoring max %d/%d\n'], ...
  context.endpointContrast.deltaNLLMaxOverSum, ...
  context.endpointContrast.nSessionsFavoringMax, ...
  numel(context.endpointContrast.sessionIndex));

end
