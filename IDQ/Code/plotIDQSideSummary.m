function plotIDQSideSummary()
% plotIDQSideSummary
%
% Two-page across-session IDQ summary, with parallel pages for the change
% and no-change sides.
%
% Reads:
%   Data/SessionAnalysis/*_sessionAnalysis.mat
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_AcrossSideSummary.mat
%   Plots/AcrossSessionSummaries/IDQSideSummary.pdf
%
% The normalized gain analysis proceeds in several stages:
%
% (1) Trial tables from all sessions are concatenated, retaining a session index so that session-specific 
% psychometric parameters can be applied to every trial. 
% (2) An initial Weibull function is fitted separately to each session. These fits provide a session-specific 
% 75%-correct coherence threshold. 
% (3) Each trial's coherence, including its noise-perturbed coherence, is divided by the threshold for its session. 
% This places all sessions on a common normalized coherence scale, where 1 corresponds to the session's estimated 
% 75%-correct threshold. 
% (4) The coherence-normalized trials from all sessions are then pooled, and a single across-session Weibull is 
% fitted to determine the common Weibull exponent and lapse rate. Thus, session-to-session differences in 
% sensitivity are represented by the session thresholds, while the shape of the normalized psychometric function is 
% estimated from all sessions together.
% (5) For each session and side, the three direction-specific noise streams are converted into trial-by-trial predictors.
% The current plotted analysis uses their step-window rectangular averages. Optional full-trial rectangular and
% leave-one-session-out kernel-weighted predictors remain available through noisePredictorType.
% 6) Only trials containing step noise are used for the gain fits. The fits are performed in the original coherence 
% units. Each session's fitted 75%-correct threshold is represented by a corresponding session-specific Weibull scale 
% parameter, calculated using the fixed across-session Weibull exponent and lapse rate. This is mathematically 
% equivalent to the earlier threshold normalization, but allows the fitted gains to remain in coherence units.
% (7) The noise predictors are added to the physical coherence step to define an effective coherence, 
% (c_eff=cstep+Xg, with negative effective coherences clipped to zero. 
% (8) The same four descriptive gain models are fitted separately for each side: one common gain for the sum of all three streams, three gains defined by 
% absolute direction, three gains defined relative to the drift direction, and separate gains for drift 
% and the summed non-drift streams. Within each model, the gain parameters are shared across sessions and are 
% estimated jointly from all eligible trials by maximizing the Bernoulli likelihood. These are marginal side-specific
% fits, not a joint six-stream model. The session-specific Weibull 
% scale parameters and the common Weibull exponent and lapse are held fixed during this optimization. 
% (9) Gains are constrained to the range −5 to 5. Approximate standard errors are obtained from the inverse of a 
% finite-difference Hessian of the negative log likelihood at the optimum, provided that the Hessian is positive 
% definite. The plotted 95% confidence intervals are the resulting Wald intervals, (g\pm1.96,SE).
%
% Noise predictors retain the sign of physical directional coherence on
% both sides. Thus, behaviorally adverse no-change-side noise is expected
% to have a negative fitted gain; no predictor or gain sign is reversed.

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
    error('plotIDQSideSummary:NoSessionAnalysisFiles', 'No sessionAnalysis files found in %s.', sessionAnalysisFolder);
end

sessionAnalyses = cell(numel(files), 1);
for iFile = 1:numel(files)
    filePath = fullfile(files(iFile).folder, files(iFile).name);
    load(filePath, 'sessionAnalysis');
    sessionAnalyses{iFile} = sessionAnalysis;
end
validateSideInputs(sessionAnalyses);

nDirs = sessionAnalyses{1}.sessionHeader.nDirs;
dirStepDeg = 360 / nDirs;
absDirLabels = strings(1, nDirs);
for d = 1:nDirs
  absDirLabels(d) = sprintf('%.0f°-%.0f°', (d - 1) * dirStepDeg, d * dirStepDeg - 1);
end
dirOffset = 360 / nDirs;
relDirLabels = {'Drift Direction', sprintf('Drift - %.0f°', dirOffset), sprintf('Drift + %.0f°', dirOffset)};

% fit the across session psychometric functions
trialTable = concatenateTrialTables(sessionAnalyses);
sessionRecords = makeSessionRecords(sessionAnalyses, 'change');

[psych, trialTable] = fitPsychometric(trialTable, [1, 2, 3], 'All Directions');
sessionRecords.threshold75 = psych.sessionFits.threshold;
sessionRecords.noisyStepAlignedCoh = psych.sessionFits.noisyStepCoh ./ psych.sessionFits.threshold;
psychometric = computeAcrossPsychometric(trialTable);

% fprintf('Session 75%% thresholds (coherence): %s\n', sprintf('%.2f ', sessionRecords.threshold75));
% fprintf(['75%% threshold across sessions: mean %.2f, SD %.2f, CV %.2f, ' ...
%   'median %.2f, range %.2f–%.2f (n = %d)\n'], ...
%   mean(sessionRecords.threshold75, 'omitnan'), ...
%   std(sessionRecords.threshold75, 'omitnan'), ...
%   std(sessionRecords.threshold75, 'omitnan') / mean(sessionRecords.threshold75, 'omitnan'), ...
%   median(sessionRecords.threshold75, 'omitnan'), ...
%   min(sessionRecords.threshold75, [], 'omitnan'), ...
%   max(sessionRecords.threshold75, [], 'omitnan'), ...
%   sum(isfinite(sessionRecords.threshold75)));

psychByDir = cell(1, nDirs);
dirTrialTables = cell(1, nDirs);
for d = 1:nDirs
  [psychByDir{d}, dirTrialTables{d}] = fitPsychometric(trialTable, d, absDirLabels(d));
end

% fit the across session regressions
% Predictor options: 'rectStep', 'rectFull', 'looStep', 'looFull'
noisePredictorType = 'rectStep';
% noisePredictorType = 'rectFull';
% noisePredictorType = 'looStep';
% noisePredictorType = 'looFull';
targetPerformance = 0.75;

sideNames = {'change', 'noChange'};
sideTitles = {'Change Side', 'No-Change Side'};
sideSummaries = struct();
for iSide = 1:numel(sideNames)
  sideName = sideNames{iSide};
  sideTitle = sideTitles{iSide};

  sideSummary = struct();
  sideSummary.sideName = sideName;
  sideSummary.sideTitle = sideTitle;
  sideSummary.kernel = computeAllDirKernel(sessionAnalyses, sideName, 'All Directions');
  sideSummary.absKernels = computeDirKernels(sessionAnalyses, nDirs, true, ...
    absDirLabels, sideName);
  sideSummary.relKernels = computeDirKernels(sessionAnalyses, nDirs, false, ...
    relDirLabels, sideName);

  sideTrialTable = selectNoisePredictor(trialTable, sessionAnalyses, ...
    noisePredictorType, sideName);
  sideSummary.noiseGain = fitIDQNoiseGain(sideTrialTable, psych.sessionFits, ...
    psych.alignedWeibull, targetPerformance);
  sideSummaries.(sideName) = sideSummary;
end

% fprintf('Noise predictor: %s\n', noisePredictorType);
% fprintf('Combined gain model NLL: %.6f\n', getFitNLL(noiseGain.combined));
% fprintf('Absolute gain model NLL: %.6f\n', getFitNLL(noiseGain.absolute));
% fprintf('Drift-relative gain model NLL: %.6f\n', getFitNLL(noiseGain.driftRelative));
% fprintf('Drift/non-drift gain model NLL: %.6f\n', getFitNLL(noiseGain.driftNonDrift));

acrossSummary = struct();
acrossSummary.noisePredictorType = noisePredictorType;
acrossSummary.sessionRecords = sessionRecords;
acrossSummary.trialTable = trialTable;
acrossSummary.createdBy = mfilename;
acrossSummary.createdAt = datetime('now');
acrossSummary.sessionAnalysisFolder = sessionAnalysisFolder;
acrossSummary.nSessions = numel(sessionAnalyses);
acrossSummary.nTrials = height(trialTable);
acrossSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);

acrossSummary.psychFit = psych;
acrossSummary.psychometric = psychometric;
acrossSummary.psychByDir = psychByDir;
acrossSummary.dirTrialTables = dirTrialTables;

acrossSummary.sessionRecords = sessionRecords;
acrossSummary.trialTable = trialTable;
acrossSummary.side = sideSummaries;

summaryMatFile = fullfile(summaryFolder, 'IDQ_AcrossSideSummary.mat');
save(summaryMatFile, 'acrossSummary', '-v7.3');

gainYLim = getNoiseGainYLim( ...
  acrossSummary.side.change.noiseGain, acrossSummary.side.noChange.noiseGain);
kernelYLim = getKernelYLim( ...
  acrossSummary.side.change, acrossSummary.side.noChange);
pdfFile = fullfile(plotFolder, 'IDQSideSummary.pdf');

fig = plotSummary(acrossSummary, 'change', gainYLim, kernelYLim, 500);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');

fig = plotSummary(acrossSummary, 'noChange', gainYLim, kernelYLim, 501);
exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'Append', true);

end

%% -------------------------------------------------------------------------
function validateSideInputs(sessionAnalyses)

requiredSAFields = { ...
  'changeSumNoiseByFrameTrial', 'changeDirNoiseByFrameTrial', 'changeSumNoiseKernel', ...
  'noChangeSumNoiseByFrameTrial', 'noChangeDirNoiseByFrameTrial', 'noChangeSumNoiseKernel'};
requiredTableVars = cell(1, 12);
iVar = 0;
for temporalPrefix = {'rectStep', 'rectFull'}
  for sideName = {'change', 'noChange'}
    names = makeRectPredictorNames(temporalPrefix{1}, sideName{1});
    for iDir = 1:3
      iVar = iVar + 1;
      requiredTableVars{iVar} = names{iDir};
    end
  end
end

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  missingSA = requiredSAFields(~isfield(SA, requiredSAFields));
  missingTable = requiredTableVars(~ismember(requiredTableVars, SA.trialTable.Properties.VariableNames));
  if ~isempty(missingSA) || ~isempty(missingTable)
    error('plotIDQSideSummary:OldSessionAnalysis', ...
      ['Session %s lacks the side-specific inventory fields. Rerun ' ...
       'makeIDQSessionAnalyses before calling plotIDQSideSummary.'], SA.fileName);
  end
end

end

%% -------------------------------------------------------------------------
function [psych, dirTrialTable] = fitPsychometric(trialTable, dirIndices, title)
% Fit selected drift directions with a psychometric function

initialBetaWeibull = 2;
initialLapse = 0.02;
targetPerformance = 0.75;
lapseBounds = [0 0.05];

dirTrialTable = trialTable(ismember(trialTable.dirIndex, dirIndices), :);
psychDir = fitIDQInitialSessionThresholds(dirTrialTable, initialBetaWeibull, initialLapse, targetPerformance);
alignedWeibull = fitIDQAcrossAlignedWeibull(psychDir.alignedCoh, dirTrialTable.correct, ...
  targetPerformance, lapseBounds);
dirTrialTable.alignedCoh = psychDir.alignedCoh;
dirTrialTable.noisyStepAlignedCoh = psychDir.noisyStepAlignedCoh;

psych = struct();
psych.title = title;
psych.initialBetaWeibull = initialBetaWeibull;
psych.initialLapse = initialLapse;
psych.targetPerformance = targetPerformance;
psych.lapseBounds = lapseBounds;
psych.sessionFits = psychDir.sessionFits;
psych.alignedWeibull = alignedWeibull;

end

%% -------------------------------------------------------------------------
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

% Put sessionIndex near the front for readability.
trialTable = movevars(trialTable, 'sessionIndex', 'Before', 1);

end

%% -------------------------------------------------------------------------
function sessionRecords = makeSessionRecords(sessionAnalyses, sideName)

nSessions = numel(sessionAnalyses);

fileName = strings(nSessions, 1);
nTrials = nan(nSessions, 1);
nStepNoiseTrials = nan(nSessions, 1);
noisyStepCoh = nan(nSessions, 1);
meanCorrect = nan(nSessions, 1);
meanCorrectStepNoise = nan(nSessions, 1);
meanRectNoise = nan(nSessions, 1);
stdRectNoise = nan(nSessions, 1);
kernelStepIntegral = nan(nSessions, 1);
kernelStepMean = nan(nSessions, 1);
kernelStepPeak = nan(nSessions, 1);

for iSession = 1:nSessions
  SA = sessionAnalyses{iSession};
  T = SA.trialTable;
  fileName(iSession) = string(SA.fileName);
  nTrials(iSession) = SA.nTrials;
  nStepNoiseTrials(iSession) = SA.nStepNoiseTrials;
  noisyStepCoh(iSession) = SA.noisyStepCoh;
  meanCorrect(iSession) = mean(T.correct, 'omitnan');
  idxNoise = T.hasStepNoise;
  meanCorrectStepNoise(iSession) = mean(T.correct(idxNoise), 'omitnan');
  fields = getSideFieldNames(sideName);
  k = SA.(fields.sumKernel).kernel;
  kStep = k(SA.stepFrames);
  kernelStepIntegral(iSession) = sum(kStep, 'omitnan');
  kernelStepMean(iSession) = mean(kStep, 'omitnan');
  kernelStepPeak(iSession) = max(abs(kStep), [], 'omitnan');
end

sessionRecords = table(fileName, nTrials, nStepNoiseTrials, noisyStepCoh, meanCorrect, ...
    meanCorrectStepNoise, meanRectNoise, stdRectNoise, kernelStepIntegral, kernelStepMean, kernelStepPeak);

end

%% -------------------------------------------------------------------------
function fig = plotSummary(acrossSummary, sideName, gainYLim, kernelYLim, figNum)

sideSummary = acrossSummary.side.(sideName);

fig = figure(figNum);
clf(fig);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 4, 4, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Summary: %s (%d sessions)', ...
    sideSummary.sideTitle, acrossSummary.nSessions), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

textAx = nexttile(tl, 12);
plotTextSummary(textAx, acrossSummary, sideName);
allColor = [0.5, 0.5, .5];
absColor = [0.0, 0.5, 0.8];
relColor = [0.8, 0.5, 1.0];

% psychometric functions

psychAx = gobjects(4,1);
for d = 1:3
  psychAx(d) = nexttile(tl, d);
  plotAlignedWeibull(psychAx(d),  acrossSummary.dirTrialTables{d}, acrossSummary.psychByDir{d});
end

psychAx(4) = nexttile(tl, 4);
plotAlignedWeibull(psychAx(4), acrossSummary.trialTable, acrossSummary.psychFit);

yl = cell2mat(get(psychAx, 'YLim'));
sharedYLim = [min(yl(:,1)), max(yl(:,2))];
set(psychAx, 'YLim', sharedYLim);

% kernels

kernelAx = gobjects(7,1);
stepPatch = gobjects(7,1);
for d = 1:3
  kernelAx(d) = nexttile(tl, d + 4);
  stepPatch(d) = plotKernel(kernelAx(d), sideSummary.absKernels{d}, ...
    absColor, false, true, sideSummary.sideTitle);
end

kernelAx(4) = nexttile(tl, 8);
stepPatch(4) = plotKernel(kernelAx(4), sideSummary.kernel, ...
  allColor, false, true, sideSummary.sideTitle);

for d = 1:3
  kernelAx(d + 4) = nexttile(tl, d+  8);
  stepPatch(d + 4) = plotKernel(kernelAx(d + 4), ...
    sideSummary.relKernels{d}, relColor, d == 1, false, sideSummary.sideTitle);
end

set(kernelAx, 'YLim', kernelYLim);
for k = 1:7
  stepPatch(k).YData = kernelYLim([1 1 2 2]);
end

% Common limits across both pages permit direct side comparisons.
plotGainBars(nexttile(tl, 13), sideSummary.noiseGain.driftRelative, 'Drift-Rel. Gains', gainYLim, relColor);
plotGainBars(nexttile(tl, 14), sideSummary.noiseGain.absolute, 'Abs. Direction Gains', gainYLim, absColor);
plotGainBars(nexttile(tl, 15), sideSummary.noiseGain.combined, 'Combined Gain', gainYLim, allColor);
plotGainBars(nexttile(tl, 16), sideSummary.noiseGain.driftNonDrift, 'Drift v. Non-Drift Gains', gainYLim, relColor);

end

%% -------------------------------------------------------------------------
function plotTextSummary(ax, acrossSummary, sideName)

sideSummary = acrossSummary.side.(sideName);

axis(ax, 'off');
txt = {
    sprintf('Created: %s', string(acrossSummary.createdAt, 'dd/MM/yyyy HH:mm'))
    sprintf('Side: %s', sideSummary.sideTitle)
    sprintf('Noise predictor: %s', acrossSummary.noisePredictorType)
    sprintf('Combined:\n  NLL %.2f, g %.3f', ...
      getFitNLL(sideSummary.noiseGain.combined), sideSummary.noiseGain.combined.gain(1))
    sprintf('Drift/non-drift:\n  NLL %.2f, g_D %.3f, g_N %.3f', ...
      getFitNLL(sideSummary.noiseGain.driftNonDrift), sideSummary.noiseGain.driftNonDrift.gain(1), ...
      sideSummary.noiseGain.driftNonDrift.gain(2))
    'Predictor signs = physical coherence signs'
    ''
    'Combined:'
    sprintf(' c_{eff} = c_{step} + g_{all}(n_1 + n_2 + n_3)');
    'Absolute:'
    sprintf(' c_{eff} = c_{step} + g_1n_1 + g_2n_2 + g_3n_3');
    'Drift-Relative:'
    sprintf(' c_{eff} = c_{step} + g_Dn_D + g_{+120}n_{+120} + g_{-120}n_{-120}');
    'Drift v. Non-Drift:'
    sprintf(' c_{eff} = c_{step} + g_Dn_D + g_N(n_{+120} + n_{-120})');
};

text(ax, -0.3, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', 'FontSize', 6.5);
end

%% -------------------------------------------------------------------------
function stepPatch = plotKernel(ax, kernel, plotColor, showLegend, showXLabel, sideTitle)

if nargin < 6
  sideTitle = 'Change Side';
end
if nargin < 5
  showXLabel = true;
end
if nargin < 4
  showLegend = false;
end
hold(ax, 'on');
plot(ax, kernel.tMS, zeros(size(kernel.tMS)), ':', 'HandleVisibility', 'off');
x1 = kernel.tMS(kernel.stepFrames(1));
x2 = kernel.tMS(kernel.stepFrames(end));

% Plot first to establish y-limits.
plot(ax, kernel.tMS, kernel.meanDiff, '-', 'LineWidth', 1.2,  'HandleVisibility', 'off');
plot(ax, kernel.tMS, kernel.rectReference, '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');
yl = ylim(ax);
stepPatch = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.7 0.7 0.7],  'EdgeColor', 'none', ...
    'FaceAlpha', 0.35, 'HandleVisibility', 'off');

% Replot on top of patch.
plot(ax, kernel.tMS, kernel.meanDiff, '-', 'Color', plotColor, 'LineWidth', 1.2, 'DisplayName', 'Kernel');
plot(ax, kernel.tMS, kernel.rectReference, '-k',  'LineWidth', 1.0, 'DisplayName', 'Mean Pre/Post Step');
if showXLabel 
  xlabel(ax, 'Trial Time (ms)');
end
ylabel(ax, sprintf('%s Kernel', sideTitle));
title(ax, kernel.title, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');

txt = sprintf([ 'PreStep Mean %.2f%%\n' 'Step Mean %.2f%%\n'], kernel.preStepMean,  kernel.stepMean);
text(ax, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontSize', 8);
xticks(ax, 0:250:1000);
xticklabels(ax, {'0', '', '', '', '1000'});
xtickangle(ax, 0);
if showLegend
  legend(ax, 'Location', 'southeast', 'FontSize', 7);
end
end

%% -------------------------------------------------------------------------
function plotAlignedWeibull(ax, T, psych)
hold(ax, 'on'); 
fit = psych.alignedWeibull;
x = T.alignedCoh; 
correct = T.correct; 

% Bin aligned coherence for plotting. 
binEdges = linspace(min(x), max(x), 10); 
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
  'MarkerFaceColor', 'k', 'MarkerSize', 4, 'DisplayName', 'binned data');
for i = find(idxPlot(:))' 
  text(ax, binCenters(i), pCorrect(i), sprintf(' %d', nTrials(i)), 'FontSize', 7, 'VerticalAlignment', 'bottom'); 
end 
xGrid = linspace(0, max(x) * 1.05, 300); 
pFit = idqWeibullP(xGrid, fit.alpha, fit.betaWeibull, fit.lapse); 
plot(ax, xGrid, pFit, '-', 'LineWidth', 1.5);
xline(ax, fit.threshold, ':', 'DisplayName', sprintf('%.0f%% threshold = %.2f', ...
  100 * fit.thresholdPerformance, fit.threshold)); 
yline(ax, fit.thresholdPerformance, ':', 'HandleVisibility', 'off'); 
xlabel(ax, 'Norm. Coh. Step'); 
ylabel(ax, 'P(hit)'); 
title(ax, psych.title); 
ylim(ax, [0.95 * min(pCorrect(idxPlot)), 1.02]);
xlim(ax, [0 max(xGrid)]); 
grid(ax, 'on'); 
box(ax, 'off'); 

normCoh = median(T.stepCoh ./ T.alignedCoh, 'omitnan');
thresholdCoh = fit.threshold * normCoh;
txt = sprintf('%.0f trials\n%.1f%% correct\n1 = %.1f%% coh.\nthresh. = %.1f%%\n\\beta=%.2f\n\\lambda=%.3f', ...
  height(T), 100 * mean(T.correct), normCoh, thresholdCoh, fit.betaWeibull, fit.lapse);
text(ax, 0.98, 0.02, txt, 'Units', 'normalized', ...
  'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
end

%% -------------------------------------------------------------------------
function plotGainBars(ax, fit, plotTitle, gainYLim, plotColor)

gains = fit.gain(:);
CI95 = fit.CI95;
xLabels = fit.predictorNames;

nGains = numel(gains);
x = 1:nGains;
bar(ax, x, gains, 'BarWidth', 0.55, 'FaceColor', plotColor, 'EdgeColor', 'k');
hold(ax, 'on');
for iGain = 1:nGains
  plot(ax, [x(iGain) x(iGain)], CI95(iGain, :), 'k-', 'LineWidth', 1.3);
end
plot(ax, x, gains, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineStyle', 'none');
yline(ax, 0, 'k:','HandleVisibility', 'off');
set(ax, 'XTick', x, 'XTickLabel', xLabels);

% A common three-position axis gives identical physical bar widths in
% the one-, two-, and three-parameter plots.
xlim(ax, [0.25 3.75]);
ylim(ax, gainYLim);
ylabel(ax, 'Noise gain');
title(ax, plotTitle, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');
txt = '';
for d = 1:numel(gains)
  txt = [txt, sprintf('%s: %.2f\n', fit.predictorNames{d}, gains(d))]; %#ok<AGROW>
end
text(ax, 0.98, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right', ...
  'FontSize', 8);
end

%% -------------------------------------------------------------------------
function gainYLim = getNoiseGainYLim(varargin)

allValues = 0;

for iSide = 1:nargin
  noiseGain = varargin{iSide};
  fits = {noiseGain.combined noiseGain.absolute noiseGain.driftRelative noiseGain.driftNonDrift};
  for iFit = 1:numel(fits)
    fit = fits{iFit};
    allValues = [allValues 
      fit.gain(:)
      fit.CI95(:)]; %#ok<AGROW>
  end
end

allValues = allValues(isfinite(allValues));
if isempty(allValues)
  gainYLim = [-1 1];
  return
end

yMin = min(allValues);
yMax = max(allValues);
yRange = yMax - yMin;
if yRange == 0
  yRange = max(1, abs(yMax));
end
padding = 0.10 * yRange;
gainYLim = [min(0, yMin - padding) yMax + padding];

end

%% -------------------------------------------------------------------------
function kernelYLim = getKernelYLim(varargin)

allValues = 0;
for iSide = 1:nargin
  sideSummary = varargin{iSide};
  kernels = [{sideSummary.kernel}, sideSummary.absKernels, sideSummary.relKernels];
  for iKernel = 1:numel(kernels)
    allValues = [allValues; kernels{iKernel}.meanDiff(:)]; %#ok<AGROW>
  end
end

allValues = allValues(isfinite(allValues));
if isempty(allValues)
  kernelYLim = [-1 1];
  return
end

yMin = min(allValues);
yMax = max(allValues);
yRange = yMax - yMin;
if yRange == 0
  yRange = max(1, abs(yMax));
end
padding = 0.10 * yRange;
kernelYLim = [min(0, yMin - padding), max(0, yMax + padding)];

end

%% -------------------------------------------------------------------------
function kernel = computeAllDirKernel(sessionAnalyses, sideName, title)

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
allCorrect = [];
allError = [];
fields = getSideFieldNames(sideName);

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  if numel(SA.tMS) ~= numel(tMS) || any(SA.tMS(:) ~= tMS(:))
    error('plotIDQSideSummary:TimeVectorMismatch', 'Session %s has a different tMS vector.', SA.fileName);
  end
  if numel(SA.stepFrames) ~= numel(stepFrames) || any(SA.stepFrames(:) ~= stepFrames(:))
    error('plotIDQSideSummary:StepFrameMismatch', 'Session %s has different stepFrames.', SA.fileName);
  end

  T = SA.trialTable;
  idxUse = T.hasStepNoise;
  correctUse = idxUse & T.correct;
  errorUse = idxUse & ~T.correct;
  allCorrect = [allCorrect, SA.(fields.sumByFrame)(:, correctUse)]; %#ok<AGROW>
  allError = [allError, SA.(fields.sumByFrame)(:, errorUse)]; %#ok<AGROW>
end

kernel = packageKernel(tMS, stepFrames, allCorrect, allError, title);

end

%% ------------------------------------------------------------------------
function kernels = computeDirKernels(sessionAnalyses, nDirs, absolute, titles, sideName)

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
fields = getSideFieldNames(sideName);
kernels = cell(1, nDirs);
nCorrect = nan(nDirs, 1);
nError = nan(nDirs, 1);

for iDir = 1:nDirs
  allCorrect = [];
  allError = [];

  for iSession = 1:numel(sessionAnalyses)
    SA = sessionAnalyses{iSession};
    T = SA.trialTable;
    idxUse = T.hasStepNoise;
    driftDirIndex = T.dirIndex(:);
    dirIndex = iDir;
    nTrial = height(T);
    sideDirNoise = SA.(fields.dirByFrame);
    noiseThisDir = nan(numel(tMS), nTrial);
    for iTrial = 1:nTrial
      if ~absolute
        if iDir == 1
          dirIndex = driftDirIndex(iTrial);
        elseif iDir == 2
          dirIndex = mod(driftDirIndex(iTrial), 3) + 1;
        else
          dirIndex = mod(driftDirIndex(iTrial) - 2, 3) + 1;
        end
      end
      noiseThisDir(:, iTrial) = sideDirNoise(:, iTrial, dirIndex);
    end
    allCorrect = [allCorrect, noiseThisDir(:, idxUse & T.correct)]; %#ok<AGROW>
    allError = [allError, noiseThisDir(:, idxUse & ~T.correct)]; %#ok<AGROW>
  end
  nCorrect(iDir) = size(allCorrect, 2);
  nError(iDir) = size(allError, 2);

  kernels{iDir} = packageKernel(tMS, stepFrames, allCorrect, allError, titles(iDir));
end

end

%%-----------------------------------------------------------------------------
function kernel = packageKernel(tMS, stepFrames, allCorrect, allError, title)

meanCorrect = mean(allCorrect, 2, 'omitnan');
meanError = mean(allError, 2, 'omitnan');

kernel = struct();
kernel.title = title;
kernel.tMS = tMS;
kernel.stepFrames = stepFrames;
kernel.meanCorrect = meanCorrect;
kernel.meanError = meanError;
kernel.meanDiff = meanCorrect - meanError;

kernel.nCorrect = size(allCorrect, 2);
kernel.nError = size(allError, 2);
kernel.nTrials = kernel.nCorrect + kernel.nError;

kernel.stepIntegral = sum(kernel.meanDiff(stepFrames), 'omitnan');
kernel.stepMean = mean(kernel.meanDiff(stepFrames), 'omitnan');
kernel.stepSD = std(kernel.meanDiff(stepFrames), 'omitnan');
kernel.stepPeak = max(abs(kernel.meanDiff(stepFrames)), [], 'omitnan');

preStepFrames = find(tMS < tMS(stepFrames(1)));
kernel.preStepFrames = preStepFrames;

kernel.preStepMean = mean(kernel.meanDiff(preStepFrames), 'omitnan');
kernel.preStepSD = std(kernel.meanDiff(preStepFrames), 'omitnan');
kernel.stepMinusPreMean = kernel.stepMean - kernel.preStepMean;

if kernel.preStepSD > 0
  kernel.stepMeanZPreSD = kernel.stepMinusPreMean / kernel.preStepSD;
else
  kernel.stepMeanZPreSD = NaN;
end
kernel.rectReference = ones(size(kernel.meanDiff)) .* kernel.preStepMean;
kernel.rectReference(stepFrames) = kernel.stepMean;

end

%% -------------------------------------------------------------------------
function psychometric = computeAcrossPsychometric(trialTable)

cohLevels = unique(trialTable.stepCoh);
cohLevels = cohLevels(:);

nLevels = numel(cohLevels);

nTrials = nan(nLevels, 1);
nCorrect = nan(nLevels, 1);
pCorrect = nan(nLevels, 1);

nTrialsNoise = nan(nLevels, 1);
nCorrectNoise = nan(nLevels, 1);
pCorrectNoise = nan(nLevels, 1);

nTrialsNoNoise = nan(nLevels, 1);
nCorrectNoNoise = nan(nLevels, 1);
pCorrectNoNoise = nan(nLevels, 1);

for iLevel = 1:nLevels
  coh = cohLevels(iLevel);

  idx = trialTable.stepCoh == coh;
  nTrials(iLevel) = sum(idx);
  nCorrect(iLevel) = sum(trialTable.correct(idx));
  pCorrect(iLevel) = mean(trialTable.correct(idx), 'omitnan');

  idxNoise = idx & trialTable.hasStepNoise;
  nTrialsNoise(iLevel) = sum(idxNoise);
  nCorrectNoise(iLevel) = sum(trialTable.correct(idxNoise));
  pCorrectNoise(iLevel) = mean(trialTable.correct(idxNoise), 'omitnan');

  idxNoNoise = idx & ~trialTable.hasStepNoise;
  nTrialsNoNoise(iLevel) = sum(idxNoNoise);
  nCorrectNoNoise(iLevel) = sum(trialTable.correct(idxNoNoise));
  pCorrectNoNoise(iLevel) = mean(trialTable.correct(idxNoNoise), 'omitnan');
end

psychometric = table( ...
  cohLevels, ...
  nTrials, nCorrect, pCorrect, ...
  nTrialsNoise, nCorrectNoise, pCorrectNoise, ...
  nTrialsNoNoise, nCorrectNoNoise, pCorrectNoNoise);

psychometric.Properties.VariableNames{1} = 'stepCoh';

end

%%------------------------------------------------------------------
function fields = getSideFieldNames(sideName)

switch sideName
  case 'change'
    fields.sumByFrame = 'changeSumNoiseByFrameTrial';
    fields.dirByFrame = 'changeDirNoiseByFrameTrial';
    fields.sumKernel = 'changeSumNoiseKernel';
    fields.predictorToken = 'Change';
  case 'noChange'
    fields.sumByFrame = 'noChangeSumNoiseByFrameTrial';
    fields.dirByFrame = 'noChangeDirNoiseByFrameTrial';
    fields.sumKernel = 'noChangeSumNoiseKernel';
    fields.predictorToken = 'NoChange';
  otherwise
    error('plotIDQSideSummary:BadSideName', 'Unknown side name: %s', sideName);
end

end

%%------------------------------------------------------------------
function sourceNames = makeRectPredictorNames(temporalPrefix, sideName)

fields = getSideFieldNames(sideName);
sourceNames = cell(1, 3);
for iDir = 1:3
  sourceNames{iDir} = sprintf('%s%sNoisePredDir%d', ...
    temporalPrefix, fields.predictorToken, iDir);
end

end

%%------------------------------------------------------------------
function trialTable = selectNoisePredictor(trialTable, sessionAnalyses, predictorType, sideName)

switch predictorType
  case 'rectStep'
    sourceNames = makeRectPredictorNames('rectStep', sideName);
    for iDir = 1:3
      trialTable.(sprintf('noisePredDir%d', iDir)) = trialTable.(sourceNames{iDir});
    end
  case 'rectFull'
    sourceNames = makeRectPredictorNames('rectFull', sideName);
    for iDir = 1:3
      trialTable.(sprintf('noisePredDir%d', iDir)) = trialTable.(sourceNames{iDir});
    end
  case {'looStep', 'looFull'}
    temporalWindow = erase(predictorType, 'loo');
    temporalWindow = lower(temporalWindow);
    looPredictor = computeIDQLOONoisePredictor(sessionAnalyses, temporalWindow, sideName);
    trialTable.noisePredDir1 = looPredictor.noisePredDir1;
    trialTable.noisePredDir2 = looPredictor.noisePredDir2;
    trialTable.noisePredDir3 = looPredictor.noisePredDir3;
  otherwise
    error('plotIDQSideSummary:BadNoisePredictorType', ...
      'Unknown noisePredictorType: %s', predictorType);
end
end

%%------------------------------------------------------------------
function nll = getFitNLL(fit)

candidateNames = {'NLL', 'nll', 'negLogLikelihood'};
for iName = 1:numel(candidateNames)
  if isfield(fit, candidateNames{iName})
    nll = fit.(candidateNames{iName});
    return
  end
end
error('plotIDQSideSummary:MissingNLL', 'Gain fit does not contain an NLL field.');
end

%%------------------------------------------------------------------
function looPredictor = computeIDQLOONoisePredictor(sessionAnalyses, temporalWindow, sideName)
% computeIDQLOONoisePredictor
%
% Compute leave-one-session-out kernel-weighted signed-noise predictors.
% A common temporal kernel is estimated from the summed noise of all other
% sessions and then applied separately to each of the three direction streams.
%
% Outputs are aligned with the concatenated trial-table order used by
% makeIDQAcrossSessionSummary:
%   noisePredDir1, noisePredDir2, noisePredDir3

nSessions = numel(sessionAnalyses);

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
nFrames = numel(tMS);
fields = getSideFieldNames(sideName);

predCell = cell(nSessions, 3);

sessionIndex = nan(nSessions, 1);
fileName = strings(nSessions, 1);
kernelStepSum = nan(nSessions, 1);
kernelStepMean = nan(nSessions, 1);
kernelStepMin = nan(nSessions, 1);
kernelStepMax = nan(nSessions, 1);
fracNegativeWeights = nan(nSessions, 1);
effectiveNFrames = nan(nSessions, 1);

weightsBySession = nan(nFrames, nSessions);

for iSession = 1:nSessions
    SA = sessionAnalyses{iSession};
    if numel(SA.tMS) ~= numel(tMS) || any(SA.tMS(:) ~= tMS(:))
        error('computeIDQLOONoisePredictor:TimeVectorMismatch', ...
            'Session %s has a different tMS vector.', SA.fileName);
    end
    if numel(SA.stepFrames) ~= numel(stepFrames) || any(SA.stepFrames(:) ~= stepFrames(:))
        error('computeIDQLOONoisePredictor:StepFrameMismatch', ...
            'Session %s has different stepFrames.', SA.fileName);
    end
    if size(SA.(fields.dirByFrame), 3) ~= 3
        error('computeIDQLOONoisePredictor:ExpectedThreeDirs', ...
            'Session %s does not contain three direction streams.', SA.fileName);
    end

    looKernel = computeLOOKernel(sessionAnalyses, iSession, sideName);
    switch temporalWindow
      case 'step'
        useFrames = stepFrames;
      case 'full'
        useFrames = (1:nFrames)';
      otherwise
        error('computeIDQLOONoisePredictor:BadTemporalWindow', ...
          'Unknown temporal window: %s', temporalWindow);
    end
    kUse = looKernel.meanDiff(useFrames);
    kSum = sum(kUse, 'omitnan');

    if ~isfinite(kSum) || kSum == 0
        error('computeIDQLOONoisePredictor:BadKernelSum', ...
            'LOO kernel sum is zero or nonfinite for session %s.', SA.fileName);
    end

    w = kUse(:) ./ kSum;

    for iDir = 1:3
        dirNoise = SA.(fields.dirByFrame)(useFrames, :, iDir);
        predCell{iSession, iDir} = sum(dirNoise .* w, 1, 'omitnan')';
    end

    sessionIndex(iSession) = iSession;
    fileName(iSession) = string(SA.fileName);
    kernelStepSum(iSession) = kSum;
    kernelStepMean(iSession) = mean(kUse, 'omitnan');
    kernelStepMin(iSession) = min(kUse, [], 'omitnan');
    kernelStepMax(iSession) = max(kUse, [], 'omitnan');
    fracNegativeWeights(iSession) = mean(w < 0);
    effectiveNFrames(iSession) = 1 / sum(w .^ 2);
    weightsBySession(useFrames, iSession) = w;
end

looPredictor = struct();
looPredictor.noisePredDir1 = vertcat(predCell{:, 1});
looPredictor.noisePredDir2 = vertcat(predCell{:, 2});
looPredictor.noisePredDir3 = vertcat(predCell{:, 3});
looPredictor.temporalWindow = temporalWindow;
looPredictor.stepFrames = stepFrames;
looPredictor.tMS = tMS;
looPredictor.weightsBySession = weightsBySession;
looPredictor.weightDiagnostics = table( ...
    sessionIndex, fileName, kernelStepSum, kernelStepMean, kernelStepMin, ...
    kernelStepMax, fracNegativeWeights, effectiveNFrames);

end

%% ------------------------------------------------------------------------
function kernel = computeLOOKernel(sessionAnalyses, leaveOutSession, sideName)

allCorrectNoise = [];
allErrorNoise = [];
fields = getSideFieldNames(sideName);

for iSession = 1:numel(sessionAnalyses)
    if iSession == leaveOutSession && numel(sessionAnalyses) > 1
        continue
    end
    SA = sessionAnalyses{iSession};
    T = SA.trialTable;

    idxUse = T.hasStepNoise;
    correctUse = idxUse & T.correct;
    errorUse = idxUse & ~T.correct;
    allCorrectNoise = [allCorrectNoise, SA.(fields.sumByFrame)(:, correctUse)]; %#ok<AGROW>
    allErrorNoise = [allErrorNoise, SA.(fields.sumByFrame)(:, errorUse)]; %#ok<AGROW>
end

kernel = struct();
kernel.meanCorrect = mean(allCorrectNoise, 2, 'omitnan');
kernel.meanError = mean(allErrorNoise, 2, 'omitnan');
kernel.meanDiff = kernel.meanCorrect - kernel.meanError;
kernel.nCorrect = size(allCorrectNoise, 2);
kernel.nError = size(allErrorNoise, 2);
kernel.nTrials = kernel.nCorrect + kernel.nError;

end
