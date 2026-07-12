function plotIDQSummary()
% plotIDQSummary
%
% Minimal first across-session IDQ summary.
%
% Reads:
%   Data/SessionAnalysis/*_sessionAnalysis.mat
%
% Writes:
%   Data/AcrossSessionSummaries/IDQ_AcrossSessionSummary.mat
%   Data/AcrossSessionSummaries/IDQ_AcrossSessionSummary.pdf
%
% This first version intentionally does not fit Weibull functions or
% regression gain models. It verifies and summarizes the reusable
% sessionAnalysis products.

domainPath = domainFolder(mfilename('fullpath'));
sessionAnalysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
summaryFolder = validFolder(fullfile(domainPath, 'Data', 'AcrossSessionSummaries'));
plotFolder = validFolder(fullfile(domainPath, 'Plots', 'AcrossSessionSummaries'));

files = dir(fullfile(sessionAnalysisFolder, '*_sessionAnalysis.mat'));
if isempty(files)
    error('plotIDQSummary:NoSessionAnalysisFiles', 'No sessionAnalysis files found in %s.', sessionAnalysisFolder);
end
fprintf('plotIDQSummary: loading %d sessionAnalysis files\n', numel(files));

sessionAnalyses = cell(numel(files), 1);
for iFile = 1:numel(files)
    filePath = fullfile(files(iFile).folder, files(iFile).name);
    load(filePath, 'sessionAnalysis');
    sessionAnalyses{iFile} = sessionAnalysis;
end

nDirs = sessionAnalyses{1}.sessionHeader.nDirs;
dirStepDeg = 360 / nDirs;
absDirLabels = strings(1, nDirs);
for d = 1:nDirs
  absDirLabels(d) = sprintf('%.0f°-%.0f°', (d - 1) * dirStepDeg, d * dirStepDeg - 1);
end
dirOffset = 360 / nDirs;
relDirLabels = {'Drift Direction', sprintf('Drift - %.0f%%', dirOffset), sprintf('Drift + %.0f%%', dirOffset)};

% compile the absolute-, relative-, and all-direction kernels
kernel = computeAllDirKernel(sessionAnalyses, 'All Directions');
absKernels = computeDirKernels(sessionAnalyses, nDirs, true, absDirLabels);
relKernels = computeDirKernels(sessionAnalyses, nDirs, false, relDirLabels);

% fit the across session psychometric functions
trialTable = concatenateTrialTables(sessionAnalyses);
sessionRecords = makeSessionRecords(sessionAnalyses);

[psych, trialTable] = fitPsychometric(trialTable, [1, 2, 3], 'All Directions');
sessionRecords.threshold75 = psych.sessionFits.threshold;
sessionRecords.noisyStepAlignedCoh = psych.sessionFits.noisyStepCoh ./ psych.sessionFits.threshold;
psychometric = computeAcrossPsychometric(trialTable);

psychByDir = cell(1, nDirs);
dirTrialTables = cell(1, nDirs);
for d = 1:nDirs
  [psychByDir{d}, dirTrialTables{d}] = fitPsychometric(trialTable, d, absDirLabels(d));
end

% rectPredictor = computeRectPredictorSummary(trialTable, sessionAnalyses);
looPredictor = computeIDQLOONoisePredictor(sessionAnalyses);
noisePredictorWeighting = 'looKernel';  % 'looKernel' or 'rectangular'
switch noisePredictorWeighting
  case 'looKernel'
    trialTable.noisePredDir1 = looPredictor.noisePredDir1;
    trialTable.noisePredDir2 = looPredictor.noisePredDir2;
    trialTable.noisePredDir3 = looPredictor.noisePredDir3;
  case 'rectangular'
    % The concatenated session trial tables already contain rectangularly
    % weighted values under the same neutral predictor names.
  otherwise
    error('makeIDQAcrossSessionSummary:BadNoisePredictorWeighting', ...
      'Unknown noisePredictorWeighting: %s', noisePredictorWeighting);
end
% trialTable.xNoiseLOO = looPredictor.xNoiseLOO;

targetPerformance = 0.75;
noiseGain = fitIDQNoiseGain(trialTable, psych.sessionFits, psych.alignedWeibull, targetPerformance);
% noiseGainLOO = fitIDQNoiseGain(trialTable, 'xNoiseLOO', psych.sessionFits, psych.alignedWeibull, ...
%   targetPerformance);
directionDiagnostics = computeIDQDirectionDiagnostics(sessionAnalyses, trialTable);
% directionDiagnostics.absolute.rectGainByDirection = fitIDQRectGainByDirection(trialTable, psych.sessionFits, ...
%   psych.alignedWeibull, targetPerformance, directionDiagnostics.absolute.directionLabels);

acrossSummary = struct();
acrossSummary.noisePredictorWeighting = noisePredictorWeighting;
acrossSummary.looPredictor = looPredictor;
acrossSummary.sessionRecords = sessionRecords;
acrossSummary.trialTable = trialTable;
acrossSummary.createdBy = mfilename;
acrossSummary.createdAt = datetime('now');
acrossSummary.sessionAnalysisFolder = sessionAnalysisFolder;
acrossSummary.nSessions = numel(sessionAnalyses);
acrossSummary.nTrials = height(trialTable);
acrossSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
acrossSummary.directionDiagnostics = directionDiagnostics;

acrossSummary.psychFit = psych;
acrossSummary.psychometric = psychometric;
acrossSummary.psychByDir = psychByDir;
acrossSummary.dirTrialTables = dirTrialTables;

acrossSummary.sessionRecords = sessionRecords;
acrossSummary.trialTable = trialTable;
acrossSummary.kernel = kernel;
acrossSummary.kernel = kernel;
acrossSummary.absKernels = absKernels;
acrossSummary.relKernels = relKernels;

% acrossSummary.rectPredictor = rectPredictor;
acrossSummary.noiseGain = noiseGain;

% acrossSummary.noiseGain.primary = 'rect';
% acrossSummary.noiseGain.rect = noiseGainRect;
% acrossSummary.noiseGain.looKernel = noiseGainLOO;

acrossSummary.directionDiagnostics = directionDiagnostics;

summaryMatFile = fullfile(summaryFolder, 'IDQ_AcrossSessionSummary.mat');
save(summaryMatFile, 'acrossSummary', '-v7.3');
fprintf('  saved %s\n', summaryMatFile);

summaryPDFFile = fullfile(plotFolder, 'IDQSessionSummary.pdf');
fig = plotSessionSummary(acrossSummary);
exportgraphics(fig, summaryPDFFile, 'ContentType', 'vector');
% fig = plotIDQDirectionDiagnosticsSummary(acrossSummary);
% exportgraphics(fig, summaryPDFFile, 'ContentType', 'vector', 'Append', true);

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
  % T.rectSumNoise = SA.rectSumNoise(:);
  % T.rectMeanNoise = SA.rectMeanNoise(:);
  tables{iSession} = T;
end
trialTable = vertcat(tables{:});

% Put sessionIndex near the front for readability.
trialTable = movevars(trialTable, 'sessionIndex', 'Before', 1);

end

%% -------------------------------------------------------------------------
function sessionRecords = makeSessionRecords(sessionAnalyses)

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
  % x = SA.rectSumNoise(idxNoise);
  % meanRectNoise(iSession) = mean(x, 'omitnan');
  % stdRectNoise(iSession) = std(x, 'omitnan');
  k = SA.sumNoiseKernel.kernel;
  kStep = k(SA.stepFrames);
  kernelStepIntegral(iSession) = sum(kStep, 'omitnan');
  kernelStepMean(iSession) = mean(kStep, 'omitnan');
  kernelStepPeak(iSession) = max(abs(kStep), [], 'omitnan');
end

sessionRecords = table(fileName, nTrials, nStepNoiseTrials, noisyStepCoh, meanCorrect, ...
    meanCorrectStepNoise, meanRectNoise, stdRectNoise, kernelStepIntegral, kernelStepMean, kernelStepPeak);

end

% %% -------------------------------------------------------------------------
% function psychometric = computeAcrossPsychometric(trialTable)
% 
% cohLevels = unique(trialTable.stepCoh);
% cohLevels = cohLevels(:);
% 
% nLevels = numel(cohLevels);
% 
% nTrials = nan(nLevels, 1);
% nCorrect = nan(nLevels, 1);
% pCorrect = nan(nLevels, 1);
% 
% nTrialsNoise = nan(nLevels, 1);
% nCorrectNoise = nan(nLevels, 1);
% pCorrectNoise = nan(nLevels, 1);
% 
% nTrialsNoNoise = nan(nLevels, 1);
% nCorrectNoNoise = nan(nLevels, 1);
% pCorrectNoNoise = nan(nLevels, 1);
% 
% for iLevel = 1:nLevels
%     coh = cohLevels(iLevel);
% 
%     idx = trialTable.stepCoh == coh;
%     nTrials(iLevel) = sum(idx);
%     nCorrect(iLevel) = sum(trialTable.correct(idx));
%     pCorrect(iLevel) = mean(trialTable.correct(idx), 'omitnan');
% 
%     idxNoise = idx & trialTable.hasStepNoise;
%     nTrialsNoise(iLevel) = sum(idxNoise);
%     nCorrectNoise(iLevel) = sum(trialTable.correct(idxNoise));
%     pCorrectNoise(iLevel) = mean(trialTable.correct(idxNoise), 'omitnan');
% 
%     idxNoNoise = idx & ~trialTable.hasStepNoise;
%     nTrialsNoNoise(iLevel) = sum(idxNoNoise);
%     nCorrectNoNoise(iLevel) = sum(trialTable.correct(idxNoNoise));
%     pCorrectNoNoise(iLevel) = mean(trialTable.correct(idxNoNoise), 'omitnan');
% end
% 
% psychometric = table(cohLevels, nTrials, nCorrect, pCorrect, nTrialsNoise, nCorrectNoise, pCorrectNoise, ...
%     nTrialsNoNoise, nCorrectNoNoise, pCorrectNoNoise);
% 
% psychometric.Properties.VariableNames{1} = 'stepCoh';
% 
% end


%% -------------------------------------------------------------------------
% function rectPredictor = computeRectPredictorSummary(trialTable, sessionAnalyses)
% 
% idxNoise = trialTable.hasStepNoise;
% 
% x = trialTable.rectSumNoise(idxNoise);
% correct = trialTable.correct(idxNoise);
% 
% rectPredictor = struct();
% 
% rectPredictor.nTrials = numel(x);
% rectPredictor.mean = mean(x, 'omitnan');
% rectPredictor.sd = std(x, 'omitnan');
% rectPredictor.meanCorrectTrials = mean(x(correct), 'omitnan');
% rectPredictor.meanErrorTrials = mean(x(~correct), 'omitnan');
% rectPredictor.meanDiffCorrectMinusError =  rectPredictor.meanCorrectTrials - rectPredictor.meanErrorTrials;
% 
% % Per-session version for diagnostics.
% nSessions = numel(sessionAnalyses);
% fileName = strings(nSessions, 1);
% nTrials = nan(nSessions, 1);
% meanX = nan(nSessions, 1);
% sdX = nan(nSessions, 1);
% meanCorrectX = nan(nSessions, 1);
% meanErrorX = nan(nSessions, 1);
% meanDiffCorrectMinusError = nan(nSessions, 1);
% 
% for iSession = 1:nSessions
%   SA = sessionAnalyses{iSession};
%   T = SA.trialTable;
%   idx = T.hasStepNoise;
% 
%   xs = SA.rectSumNoise(idx);
%   cs = T.correct(idx);
% 
%   fileName(iSession) = string(SA.fileName);
%   nTrials(iSession) = numel(xs);
%   meanX(iSession) = mean(xs, 'omitnan');
%   sdX(iSession) = std(xs, 'omitnan');
%   meanCorrectX(iSession) = mean(xs(cs), 'omitnan');
%   meanErrorX(iSession) = mean(xs(~cs), 'omitnan');
%   meanDiffCorrectMinusError(iSession) = meanCorrectX(iSession) - meanErrorX(iSession);
% end
% rectPredictor.bySession = table(fileName, nTrials, meanX, sdX, meanCorrectX, meanErrorX,  meanDiffCorrectMinusError);
% end

%% -------------------------------------------------------------------------
function fig = plotSessionSummary(acrossSummary)

fig = figure(500);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 4, 4, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Summary (%d sessions)', acrossSummary.nSessions), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

textAx = nexttile(tl, 12);
plotTextSummary(textAx, acrossSummary);

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
  stepPatch(d) = plotKernel(kernelAx(d), acrossSummary.absKernels{d});
end

kernelAx(4) = nexttile(tl, 8);
stepPatch(4) = plotKernel(kernelAx(4), acrossSummary.kernel);

for d = 1:3
  kernelAx(d + 4) = nexttile(tl, d+  8);
  stepPatch(d + 4) = plotKernel(kernelAx(d + 4), acrossSummary.relKernels{d}, d == 1, false);
end

yl = cell2mat(get(kernelAx, 'YLim'));
sharedYLim = [min(yl(:,1)), max(yl(:,2))];
set(kernelAx, 'YLim', sharedYLim);
for k = 1:7
  stepPatch(k).YData = sharedYLim([1 1 2 2]);
end


% Use a common y-axis scale for all gain plots.
gainYLim = getNoiseGainYLim(acrossSummary.noiseGain);

axGain = nexttile(tl, 15);
fit = acrossSummary.noiseGain.combined;
plotGainBars(axGain, fit.gain, fit.CI95, {'Combined'}, 'Combined three-direction gain', gainYLim);

axAbsoluteGain = nexttile(tl, 14);
fit = acrossSummary.noiseGain.absolute;
plotGainBars(axAbsoluteGain, fit.gain, fit.CI95, {'Dir 1', 'Dir 2', 'Dir 3'}, ...
  'Absolute-direction gains', gainYLim);

axRelativeGain = nexttile(tl, 13);
fit = acrossSummary.noiseGain.driftRelative;
plotGainBars(axRelativeGain, fit.gain, fit.CI95, ...
  {'Drift', '+120', '-120'}, 'Drift-relative gains', gainYLim);

axDriftNonDriftGain = nexttile(tl, 16);
fit = acrossSummary.noiseGain.driftNonDrift;
plotGainBars(axDriftNonDriftGain, fit.gain, fit.CI95, ...
  {'Drift', 'Non-drift'}, 'Drift and pooled non-drift gains', gainYLim);

end

%% -------------------------------------------------------------------------
function plotTextSummary(ax, acrossSummary)

axis(ax, 'off');

R = acrossSummary.sessionRecords;
fit = acrossSummary.psychFit.alignedWeibull;
gainFit = acrossSummary.noiseGain.combined;
txt = {
    sprintf('Plot Created: %s', string(acrossSummary.createdAt))
    sprintf('%ld Sessions, %ld Trials', acrossSummary.nSessions, acrossSummary.nTrials)
    sprintf('Mean percent correct %.1f%%', mean(R.meanCorrect, 'omitnan') * 100.0)
    ''
    sprintf('%ld Noise Trials', acrossSummary.kernel.nTrials)
    sprintf('Mean Step Coherence: %.0f%%', mean(R.noisyStepCoh, 'omitnan'))
    sprintf('%.1f%% correct (%d hit, %d miss)', 100.0 * acrossSummary.kernel.nCorrect / acrossSummary.kernel.nTrials, ...
            acrossSummary.kernel.nCorrect, acrossSummary.kernel.nError)
    ''
    sprintf('Beta & lapse: %.3f, %.3f', fit.betaWeibull, fit.lapse)
    sprintf('Normalization Threshold: %.0f%%', 100 * fit.thresholdPerformance)
    ''
    sprintf('Int. noise gain: %.2f', gainFit.gain)
    sprintf('Int. noise 95%% CI: [%.2f %.2f]', gainFit.CI95(1), gainFit.CI95(2))
    ''
    sprintf('Step mean (SD): %.2f%% (%.2f%%)', acrossSummary.kernel.stepMean, acrossSummary.kernel.stepSD);
    sprintf('PreStep mean (SD): %.2f%% (%.2f%%)', acrossSummary.kernel.preStepMean, acrossSummary.kernel.preStepSD)
    sprintf('Step z: %.2f', acrossSummary.kernel.stepMeanZPreSD)
    };

text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', 'FontSize', 6);

% sprintf('Int. noise mean (SD): %.3g (%.3g)', acrossSummary.rectPredictor.mean, acrossSummary.rectPredictor.sd)
% ''

end


%% -------------------------------------------------------------------------
function stepPatch = plotKernel(ax, kernel, showLegend, showXLabel)

if nargin < 4
  showXLabel = true;
end
if nargin < 3
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
plot(ax, kernel.tMS, kernel.meanDiff, '-b', 'LineWidth', 1.2, 'DisplayName', 'Kernel');
plot(ax, kernel.tMS, kernel.rectReference, '-k',  'LineWidth', 1.0, 'DisplayName', 'Mean Pre/Post Step');
if showXLabel 
  xlabel(ax, 'Trial Time (ms)');
end
ylabel(ax, 'Change Side Kernel');
title(ax, kernel.title, 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');

txt = sprintf([ 'PreStep Mean %.2f%%\n' 'Step Mean %.2f%%\n'], kernel.preStepMean,  kernel.stepMean);
text(ax, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontSize', 8);
xticks(ax, 0:250:1000);
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
plot(ax, xGrid, pFit, '-', 'LineWidth', 1.5, 'DisplayName', ...
  sprintf('Weibull \\beta=%.2f, \\lambda=%.3f', fit.betaWeibull, fit.lapse));
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

txt = sprintf('%.0f trials\n%.1f%% correct\n\\beta=%.2f\n\\lambda=%.3f', height(T), ...
  100 * mean(T.correct), fit.betaWeibull, fit.lapse);
text(ax, 0.98, 0.02, txt, 'Units', 'normalized', ...
  'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8);
end

%% -------------------------------------------------------------------------
% function plotSessionPerformance(ax, sessionRecords)
% 
% plot(ax, sessionRecords.meanCorrectStepNoise, 'ko', 'LineWidth', 1.2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
% xlabel(ax, 'Session');
% ylabel(ax, 'P(hit), noise step');
% title(ax, 'Noise Step Performance by Session');
% xlim(ax, [0, numel(sessionRecords.meanCorrectStepNoise) + 1]);
% ylim(ax, [0.5 1.0]);
% grid(ax, 'on');
% box(ax, 'off');
% 
% end

%% -------------------------------------------------------------------------
% function plotSessionKernelMeans(ax, sessionRecords)
% 
% hold(ax, 'on');
% plot(ax, sessionRecords.kernelStepMean, 'ko', 'LineWidth', 1.2, 'MarkerSize', 6, 'MarkerFaceColor', 'k');
% yline(ax, 0, 'k:');
% xlabel(ax, 'Session');
% ylabel(ax, 'Kernel Mean Over Step');
% ytickformat(ax, 'percentage');
% xlim(ax, [0, numel(sessionRecords.kernelStepMean) + 1]);
% title(ax, 'Session Kernel Means');
% grid(ax, 'on');
% box(ax, 'off');
% 
% end

%% -------------------------------------------------------------------------
function fig = plotIDQDirectionDiagnosticsSummary(acrossSummary)
% plotIDQDirectionDiagnosticsSummary
%
% Diagnostic page for direction-specific IDQ behavior and kernels.
% These diagnostics are not the primary readout analysis.

Dabs = acrossSummary.directionDiagnostics.absolute;
Daln = acrossSummary.directionDiagnostics.aligned;

fig = figure(501);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 11 8.5], ...
              'PaperOrientation', 'landscape', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title(tl, 'IDQ direction diagnostics', ...
  'Interpreter', 'none', ...
  'FontWeight', 'bold');

% ax = nexttile(tl, 1);
% plotDirectionBehavior(ax, Dabs.behavior);

% ax = nexttile(tl, 2);
% plotDirectionKernelSummary(ax, Dabs.kernels, 'Absolute direction kernel step means');

ax = nexttile(tl, 3);
plotRectGainByDirection(ax, Dabs.rectGainByDirection);

% ax = nexttile(tl, 4);
% plotDirectionKernels(ax, Dabs.kernels, 'Absolute changed-side kernels');
% 
% ax = nexttile(tl, 5);
% plotDirectionKernels(ax, Daln.kernels, 'Aligned changed-side kernels');

ax = nexttile(tl, 6);
plotDirectionDiagnosticsText(ax, acrossSummary);

end

%% -------------------------------------------------------------------------
function plotDirectionBehavior(ax, behavior)

bar(ax, behavior.pCorrect);
hold(ax, 'on');
yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');
set(ax, 'XTick', 1:height(behavior), 'XTickLabel', behavior.directionLabel);
xtickangle(ax, 30);
ylim(ax, [0.45 1.0]);
ylabel(ax, 'P(correct)');
title(ax, 'Behavior by physical drift direction', 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');
for i = 1:height(behavior)
    text(ax, i, behavior.pCorrect(i), sprintf('n=%d', behavior.nTrials(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
end
end

%% -------------------------------------------------------------------------
function plotDirectionKernelSummary(ax, kernels, plotTitle)

S = kernels.summaryTable;
bar(ax, S.stepMean);
hold(ax, 'on');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
set(ax, 'XTick', 1:height(S), 'XTickLabel', S.directionLabel);
xtickangle(ax, 30);
ylabel(ax, 'Step mean kernel');
title(ax, plotTitle, 'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');
for i = 1:height(S)
  y = S.stepMean(i);
  if y >= 0
    va = 'bottom';
  else
    va = 'top';
  end
  text(ax, i, y, sprintf('%.3g', y), 'HorizontalAlignment', 'center', 'VerticalAlignment', va, 'FontSize', 7);
end
yl = ylim(ax);
m = max(abs(yl));
ylim(ax, [-m m]);
end

%% -------------------------------------------------------------------------
% function plotDirectionKernels(ax, kernels, plotTitle)
% 
% hold(ax, 'on');
% tMS = kernels.tMS;
% stepFrames = kernels.stepFrames;
% 
% plot(ax, tMS, zeros(size(tMS)), 'k:', 'HandleVisibility', 'off');
% 
% x1 = tMS(stepFrames(1));
% x2 = tMS(stepFrames(end));
% 
% % Initial plot to set y limits.
% for iDir = 1:numel(kernels.summaryTable.directionLabel)
%     plot(ax, tMS, kernels.meanDiff(:, iDir), ...
%         'LineWidth', 1.0, 'DisplayName', char(kernels.summaryTable.directionLabel(iDir)), 'HandleVisibility', 'off');
% end
% 
% yl = ylim(ax);
% patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.7 0.7 0.7], ...
%         'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
% 
% % Replot on top of patch.
% for iDir = 1:numel(kernels.summaryTable.directionLabel)
%     plot(ax, tMS, kernels.meanDiff(:, iDir), ...
%         'LineWidth', 1.0, 'DisplayName', char(kernels.summaryTable.directionLabel(iDir)));
% end
% 
% xlabel(ax, 'Time from trial start (ms)');
% ylabel(ax, 'Correct - error noise');
% title(ax, plotTitle, 'Interpreter', 'none');
% 
% grid(ax, 'on');
% box(ax, 'off');
% legend(ax, 'Location', 'northwest', 'FontSize', 7);
% 
% end

%% -------------------------------------------------------------------------
function plotDirectionDiagnosticsText(ax, acrossSummary)

axis(ax, 'off');

Dabs = acrossSummary.directionDiagnostics.absolute;
Daln = acrossSummary.directionDiagnostics.aligned;

B = Dabs.behavior;
Sabs = Dabs.kernels.summaryTable;
Saln = Daln.kernels.summaryTable;

[~, bestDir] = max(B.pCorrect);
[~, worstDir] = min(B.pCorrect);

behaviorRange = max(B.pCorrect) - min(B.pCorrect);

txt = {
    sprintf('Best physical direction: %s, %.1f%%', B.directionLabel(bestDir), B.pCorrect(bestDir) * 100.0)
    sprintf('Worst physical direction: %s, %.1f%%', B.directionLabel(worstDir), B.pCorrect(worstDir) * 100.0)
    sprintf('Behavior range: %.3f', behaviorRange)
    ''
    'Absolute kernel step means:'
    sprintf('  %s: %.2g', Sabs.directionLabel(1), Sabs.stepMean(1))
    sprintf('  %s: %.2g', Sabs.directionLabel(2), Sabs.stepMean(2))
    sprintf('  %s: %.2g', Sabs.directionLabel(3), Sabs.stepMean(3))
    ''
    'Aligned kernel step means:'
    sprintf('  %s: %.2g', Saln.directionLabel(1), Saln.stepMean(1))
    sprintf('  %s: %.2g', Saln.directionLabel(2), Saln.stepMean(2))
    sprintf('  %s: %.2g', Saln.directionLabel(3), Saln.stepMean(3))
    ''
    'Diagnostic only; primary gain remains collapsed across directions.'
    };

text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', 'FontSize', 8);

end

%% -------------------------------------------------------------------------
function plotGainBars(ax, gains, CI95, xLabels, plotTitle, gainYLim)

gains = gains(:);
nGains = numel(gains);
x = 1:nGains;
bar(ax, x, gains, 'BarWidth', 0.55);
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

end

%% -------------------------------------------------------------------------
function gainYLim = getNoiseGainYLim(noiseGain)

fits = {noiseGain.combined noiseGain.absolute noiseGain.driftRelative noiseGain.driftNonDrift};

allValues = 0;

for iFit = 1:numel(fits)
  fit = fits{iFit};
  allValues = [allValues 
    fit.gain(:)
    fit.CI95(:)]; %#ok<AGROW>
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
function plotRectGainByDirection(ax, rectGainByDirection)

S = rectGainByDirection.summaryTable;
hold(ax, 'on');
x = 1:height(S);
bar(ax, x, S.gain);

for i = 1:height(S)
  plot(ax, [x(i) x(i)], [S.CI95Low(i) S.CI95High(i)], 'k-', 'LineWidth', 1.3);
end

plot(ax, x, S.gain, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
yline(ax, rectGainByDirection.flatPrediction, 'k--', 'DisplayName', 'flat prediction');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
set(ax, 'XTick', x, 'XTickLabel', S.directionLabel);
xtickangle(ax, 30);
ylabel(ax, 'Noise Gain');
title(ax, 'Gain by Drift Direction', 'Interpreter', 'none');
grid(ax, 'on');
box(ax, 'off');

yl = ylim(ax);
ylim(ax, [min([yl(1), min(S.CI95Low)-0.1, -0.1]), max([yl(2), rectGainByDirection.flatPrediction+0.1])]);

for i = 1:height(S)
  text(ax, x(i), S.CI95High(i), sprintf('zF %.1f', S.zVsFlat(i)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
end

end

%% -------------------------------------------------------------------------
function kernel = computeAllDirKernel(sessionAnalyses, title)

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
allCorrect = [];
allError = [];

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  if numel(SA.tMS) ~= numel(tMS) || any(SA.tMS(:) ~= tMS(:))
    error('plotIDQSummary:TimeVectorMismatch', 'Session %s has a different tMS vector.', SA.fileName);
  end
  if numel(SA.stepFrames) ~= numel(stepFrames) || any(SA.stepFrames(:) ~= stepFrames(:))
    error('plotIDQSummary:StepFrameMismatch', 'Session %s has different stepFrames.', SA.fileName);
  end

  T = SA.trialTable;
  idxUse = T.hasStepNoise;
  correctUse = idxUse & T.correct;
  errorUse = idxUse & ~T.correct;
  allCorrect = [allCorrect, SA.sumNoiseByFrameTrial(:, correctUse)]; %#ok<AGROW>
  allError = [allError, SA.sumNoiseByFrameTrial(:, errorUse)]; %#ok<AGROW>
end

kernel = packageKernel(tMS, stepFrames, allCorrect, allError, title);

end

%% ------------------------------------------------------------------------
function kernels = computeDirKernels(sessionAnalyses, nDirs, absolute, titles)

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;
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
    sideIndex = T.sideIndex(:)';
    driftDirIndex = T.dirIndex(:);
    dirIndex = iDir;
    nTrial = height(T);
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
      noiseThisDir(:, iTrial) = squeeze(SA.noiseBySideDir(sideIndex(iTrial), dirIndex, :, iTrial));
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