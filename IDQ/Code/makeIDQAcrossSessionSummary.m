function acrossSummary = makeIDQAcrossSessionSummary()
% makeIDQAcrossSessionSummary
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
    error('makeIDQAcrossSessionSummary:NoSessionAnalysisFiles', ...
        'No sessionAnalysis files found in %s.', sessionAnalysisFolder);
end

fprintf('makeIDQAcrossSessionSummary: loading %d sessionAnalysis files\n', numel(files));

sessionAnalyses = cell(numel(files), 1);
for iFile = 1:numel(files)
    filePath = fullfile(files(iFile).folder, files(iFile).name);
    load(filePath, 'sessionAnalysis');
    sessionAnalyses{iFile} = sessionAnalysis;
end

trialTable = concatenateTrialTables(sessionAnalyses);
sessionRecords = makeSessionRecords(sessionAnalyses);
initialBetaWeibull = 3;
initialLapse = 0.02;
targetPerformance = 0.75;
lapseBounds = [0 0.05];

% Pass 1: fixed guessed beta/lapse for session thresholds.
psychPass1 = fitIDQInitialSessionThresholds( ...
    trialTable, ...
    initialBetaWeibull, ...
    initialLapse, ...
    targetPerformance);

trialTable.alignedCoh_pass1 = psychPass1.alignedCoh;
trialTable.noisyStepAlignedCoh_pass1 = psychPass1.noisyStepAlignedCoh;

alignedWeibullPass1 = fitIDQAcrossAlignedWeibull( ...
    trialTable.alignedCoh_pass1, ...
    trialTable.correct, ...
    targetPerformance, ...
    lapseBounds);

% Pass 2: refit session thresholds using pass-1 beta/lapse.
psychPass2 = fitIDQInitialSessionThresholds( ...
    trialTable, ...
    alignedWeibullPass1.betaWeibull, ...
    alignedWeibullPass1.lapse, ...
    targetPerformance);

trialTable.alignedCoh_pass2 = psychPass2.alignedCoh;
trialTable.noisyStepAlignedCoh_pass2 = psychPass2.noisyStepAlignedCoh;

alignedWeibullPass2 = fitIDQAcrossAlignedWeibull( ...
    trialTable.alignedCoh_pass2, ...
    trialTable.correct, ...
    targetPerformance, ...
    lapseBounds);

sessionRecords.threshold75_pass1 = psychPass1.sessionFits.threshold;
sessionRecords.threshold75_pass2 = psychPass2.sessionFits.threshold;
sessionRecords.noisyStepAlignedCoh_pass1 = ...
  psychPass1.sessionFits.noisyStepCoh ./ psychPass1.sessionFits.threshold;
sessionRecords.noisyStepAlignedCoh_pass2 = ...
  psychPass2.sessionFits.noisyStepCoh ./ psychPass2.sessionFits.threshold;
sessionRecords.threshold75_changeFrac = ...
  (sessionRecords.threshold75_pass2 - sessionRecords.threshold75_pass1) ./ ...
  sessionRecords.threshold75_pass1;

psych = struct();
psych.initialBetaWeibull = initialBetaWeibull;
psych.initialLapse = initialLapse;
psych.targetPerformance = targetPerformance;
psych.lapseBounds = lapseBounds;

psych.pass1.sessionFits = psychPass1.sessionFits;
psych.pass1.alignedWeibull = alignedWeibullPass1;

psych.pass2.sessionFits = psychPass2.sessionFits;
psych.pass2.alignedWeibull = alignedWeibullPass2;

psych.finalSessionFits = psychPass2.sessionFits;
psych.finalAlignedWeibull = alignedWeibullPass2;

psychometric = computeAcrossPsychometric(trialTable);
kernel = computeAcrossSignedNoiseKernel(sessionAnalyses);
rectPredictor = computeRectPredictorSummary(trialTable, sessionAnalyses);
looPredictor = computeIDQLOONoisePredictor(sessionAnalyses);
trialTable.xNoiseLOO = looPredictor.xNoiseLOO;


noiseGainRect = fitIDQNoiseGain(trialTable, 'rectSumNoise', psych.pass2.sessionFits, ...
  psych.pass2.alignedWeibull, targetPerformance);
noiseGainLOO = fitIDQNoiseGain(trialTable, 'xNoiseLOO', psych.pass2.sessionFits, psych.pass2.alignedWeibull, ...
  targetPerformance);
directionDiagnostics = computeIDQDirectionDiagnostics(sessionAnalyses, trialTable);
directionDiagnostics.absolute.rectGainByDirection = fitIDQRectGainByDirection(trialTable, psych.pass2.sessionFits, ...
  psych.pass2.alignedWeibull, targetPerformance, directionDiagnostics.absolute.directionLabels);

acrossSummary = struct();
acrossSummary.looPredictor = looPredictor;
acrossSummary.psychFit = psych;
acrossSummary.psychometric = psychometric;
acrossSummary.noiseGain.primary = 'rect';
acrossSummary.noiseGain.rect = noiseGainRect;
acrossSummary.noiseGain.looKernel = noiseGainLOO;
acrossSummary.createdBy = mfilename;
acrossSummary.createdAt = datetime('now');
acrossSummary.sessionAnalysisFolder = sessionAnalysisFolder;
acrossSummary.nSessions = numel(sessionAnalyses);
acrossSummary.nTrials = height(trialTable);
acrossSummary.nStepNoiseTrials = sum(trialTable.hasStepNoise);
acrossSummary.directionDiagnostics = directionDiagnostics;

acrossSummary.sessionRecords = sessionRecords;
acrossSummary.trialTable = trialTable;
acrossSummary.kernel = kernel;
acrossSummary.rectPredictor = rectPredictor;

summaryMatFile = fullfile(summaryFolder, 'IDQ_AcrossSessionSummary.mat');
save(summaryMatFile, 'acrossSummary', '-v7.3');
fprintf('  saved %s\n', summaryMatFile);

summaryPDFFile = fullfile(plotFolder, 'IDQ_AcrossSessionSummary.pdf');
fig = plotIDQAcrossSessionSummary(acrossSummary);
exportgraphics(fig, summaryPDFFile, 'ContentType', 'vector');
fig = plotIDQDirectionDiagnosticsSummary(acrossSummary);
exportgraphics(fig, summaryPDFFile, 'ContentType', 'vector', 'Append', true);

end

%% -------------------------------------------------------------------------
function trialTable = concatenateTrialTables(sessionAnalyses)

tables = cell(numel(sessionAnalyses), 1);

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};
  T = SA.trialTable;

  T.sessionIndex = repmat(iSession, height(T), 1);
  T.noisyStepCoh = repmat(SA.noisyStepCoh, height(T), 1);
  T.rectSumNoise = SA.rectSumNoise(:);
  T.rectMeanNoise = SA.rectMeanNoise(:);
  % T.rectDriftMinusNonNoise = SA.rectDriftMinusNonNoise(:);

  % Compatibility alias used by older fitting/plotting code.
  % Primary IDQ flat-readout predictor is now rectSumNoise.
  T.rectNoisePredictor = T.rectSumNoise;
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

  x = SA.rectSumNoise(idxNoise);
  meanRectNoise(iSession) = mean(x, 'omitnan');
  stdRectNoise(iSession) = std(x, 'omitnan');

  k = SA.sumNoiseKernel.kernel;
  kStep = k(SA.stepFrames);

  kernelStepIntegral(iSession) = sum(kStep, 'omitnan');
  kernelStepMean(iSession) = mean(kStep, 'omitnan');
  kernelStepPeak(iSession) = max(abs(kStep), [], 'omitnan');
end

sessionRecords = table( ...
    fileName, ...
    nTrials, ...
    nStepNoiseTrials, ...
    noisyStepCoh, ...
    meanCorrect, ...
    meanCorrectStepNoise, ...
    meanRectNoise, ...
    stdRectNoise, ...
    kernelStepIntegral, ...
    kernelStepMean, ...
    kernelStepPeak);

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

%% -------------------------------------------------------------------------
function kernel = computeAcrossSignedNoiseKernel(sessionAnalyses)
% This assumes all sessions have the same tMS and stepFrames. If they do not,
% the code should crash or produce an obvious mismatch.

tMS = sessionAnalyses{1}.tMS;
stepFrames = sessionAnalyses{1}.stepFrames;

allCorrectNoise = [];
allErrorNoise = [];

for iSession = 1:numel(sessionAnalyses)
  SA = sessionAnalyses{iSession};

  if numel(SA.tMS) ~= numel(tMS) || any(SA.tMS(:) ~= tMS(:))
    error('makeIDQAcrossSessionSummary:TimeVectorMismatch', ...
      'Session %s has a different tMS vector.', SA.fileName);
  end

  if numel(SA.stepFrames) ~= numel(stepFrames) || any(SA.stepFrames(:) ~= stepFrames(:))
    error('makeIDQAcrossSessionSummary:StepFrameMismatch', ...
      'Session %s has different stepFrames.', SA.fileName);
  end

  T = SA.trialTable;
  idxUse = T.hasStepNoise;

  correctUse = idxUse & T.correct;
  errorUse = idxUse & ~T.correct;

  allCorrectNoise = [allCorrectNoise, SA.sumNoiseByFrameTrial(:, correctUse)]; %#ok<AGROW>
  allErrorNoise = [allErrorNoise, SA.sumNoiseByFrameTrial(:, errorUse)]; %#ok<AGROW>
end

meanCorrect = mean(allCorrectNoise, 2, 'omitnan');
meanError = mean(allErrorNoise, 2, 'omitnan');

kernel = struct();

kernel.tMS = tMS;
kernel.stepFrames = stepFrames;

kernel.meanCorrect = meanCorrect;
kernel.meanError = meanError;
kernel.meanDiff = meanCorrect - meanError;

kernel.nCorrect = size(allCorrectNoise, 2);
kernel.nError = size(allErrorNoise, 2);
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

kernel.rectReference = zeros(size(kernel.meanDiff));
kernel.rectReference(stepFrames) = kernel.stepMean;
end

%% -------------------------------------------------------------------------
function rectPredictor = computeRectPredictorSummary(trialTable, sessionAnalyses)

idxNoise = trialTable.hasStepNoise;

x = trialTable.rectSumNoise(idxNoise);
correct = trialTable.correct(idxNoise);

rectPredictor = struct();

rectPredictor.nTrials = numel(x);
rectPredictor.mean = mean(x, 'omitnan');
rectPredictor.sd = std(x, 'omitnan');
rectPredictor.meanCorrectTrials = mean(x(correct), 'omitnan');
rectPredictor.meanErrorTrials = mean(x(~correct), 'omitnan');
rectPredictor.meanDiffCorrectMinusError =  rectPredictor.meanCorrectTrials - rectPredictor.meanErrorTrials;

% Per-session version for diagnostics.
nSessions = numel(sessionAnalyses);
fileName = strings(nSessions, 1);
nTrials = nan(nSessions, 1);
meanX = nan(nSessions, 1);
sdX = nan(nSessions, 1);
meanCorrectX = nan(nSessions, 1);
meanErrorX = nan(nSessions, 1);
meanDiffCorrectMinusError = nan(nSessions, 1);

for iSession = 1:nSessions
  SA = sessionAnalyses{iSession};
  T = SA.trialTable;
  idx = T.hasStepNoise;

  xs = SA.rectSumNoise(idx);
  cs = T.correct(idx);

  fileName(iSession) = string(SA.fileName);
  nTrials(iSession) = numel(xs);
  meanX(iSession) = mean(xs, 'omitnan');
  sdX(iSession) = std(xs, 'omitnan');
  meanCorrectX(iSession) = mean(xs(cs), 'omitnan');
  meanErrorX(iSession) = mean(xs(~cs), 'omitnan');
  meanDiffCorrectMinusError(iSession) = meanCorrectX(iSession) - meanErrorX(iSession);
end
rectPredictor.bySession = table(fileName, nTrials, meanX, sdX, meanCorrectX, meanErrorX,  meanDiffCorrectMinusError);
end

%% -------------------------------------------------------------------------
function fig = plotIDQAcrossSessionSummary(acrossSummary)

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 14 8], 'WindowStyle', 'docked');

tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title(tl, sprintf('IDQ Across-session summary (%d sessions)', acrossSummary.nSessions), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

axText = nexttile(tl, 1);
plotAcrossTextSummary(axText, acrossSummary);

axPsych = nexttile(tl, 2);
plotAlignedPsychometricWithWeibull(axPsych, acrossSummary);

axKernel = nexttile(tl, 3);
plotAcrossKernel(axKernel, acrossSummary.kernel);

axSessPerf = nexttile(tl, 4);
plotSessionPerformance(axSessPerf, acrossSummary.sessionRecords);

axGain = nexttile(tl, 5);
plotNoiseGainSummary(axGain, acrossSummary.noiseGain);

axKernelSess = nexttile(tl, 6);
plotSessionKernelMeans(axKernelSess, acrossSummary.sessionRecords);

end

%% -------------------------------------------------------------------------
function plotAcrossTextSummary(ax, acrossSummary)

axis(ax, 'off');

R = acrossSummary.sessionRecords;
fit1 = acrossSummary.psychFit.pass1.alignedWeibull;
fit = acrossSummary.psychFit.pass2.alignedWeibull;
gainFit = acrossSummary.noiseGain.rect;
looDiag = acrossSummary.looPredictor.weightDiagnostics;
txt = {
    sprintf('Plot Created: %s', string(acrossSummary.createdAt))
    sprintf('%ld Sessions, %ld Trials', acrossSummary.nSessions, acrossSummary.nTrials)
    ''
    sprintf('Mean performance: %.1f%%', mean(R.meanCorrect, 'omitnan') * 100.0)
    sprintf('Mean noise step performance: %.1f%%', mean(R.meanCorrectStepNoise, 'omitnan') * 100.0)
    sprintf('Step coh: %.0f%%', mean(R.noisyStepCoh, 'omitnan'))
    ''
    sprintf('Rect summed-noise mean (SD): %.3g (%.3g)', acrossSummary.rectPredictor.mean, acrossSummary.rectPredictor.sd)
    ''
    sprintf('Pass1 beta/lapse: %.3f / %.3f', fit1.betaWeibull, fit1.lapse)
    sprintf('Pass2 beta/lapse: %.3f / %.3f', fit.betaWeibull, fit.lapse)
    sprintf('Normalization Threshold: %.0f%%', 100 * fit.thresholdPerformance)
    ''
    sprintf('Rect noise gain: %.3f; 95%% CI: [%.3f %.3f]', gainFit.gain, gainFit.CI95(1), gainFit.CI95(2))
    sprintf('Flat prediction: %.1f', gainFit.flatPrediction)
    sprintf('z vs flat: %.2f', gainFit.zVsFlat)
    ''
    sprintf('Across kernel step mean (SD): %.2f%% (%.2f%%)', ...
            acrossSummary.kernel.stepMean, acrossSummary.kernel.stepSD);
    sprintf('Across kernel preStep mean (SD): %.2f%% (%.2f%%)', ...
            acrossSummary.kernel.preStepMean, acrossSummary.kernel.preStepSD)
    sprintf('Across kernel step z: %.2f', acrossSummary.kernel.stepMeanZPreSD)
    ''
    sprintf('LOO Min Effective Step Frames: %.2f', min(looDiag.effectiveNFrames))
    sprintf('LOO Max Frac Neg Weights: %.2f', max(looDiag.fracNegativeWeights))
    };


text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', 'FontSize', 8);

end

%% -------------------------------------------------------------------------%% -------------------------------------------------------------------------
function plotAcrossKernel(ax, kernel)

hold(ax, 'on');

plot(ax, kernel.tMS, zeros(size(kernel.tMS)), ':', 'HandleVisibility', 'off');

x1 = kernel.tMS(kernel.stepFrames(1));
x2 = kernel.tMS(kernel.stepFrames(end));

% Plot first to establish y-limits.
plot(ax, kernel.tMS, kernel.meanDiff, '-', 'LineWidth', 1.2,  'HandleVisibility', 'off');
plot(ax, kernel.tMS, kernel.rectReference, '--', 'LineWidth', 1.0, 'HandleVisibility', 'off');

yl = ylim(ax);
patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8],  'EdgeColor', 'none', ...
    'FaceAlpha', 0.35, 'HandleVisibility', 'off');

% Replot on top of patch.
plot(ax, kernel.tMS, kernel.meanDiff, '-b', 'LineWidth', 1.2, 'DisplayName', 'correct - error');
plot(ax, kernel.tMS, kernel.rectReference, '-k',  'LineWidth', 1.0, 'DisplayName', 'step-mean rectangle');

xlabel(ax, 'Time from trial start (ms)');
ylabel(ax, 'Summed changed-side noise, correct - error');
title(ax, 'Across summed-noise kernel', 'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');

txt = sprintf(['correct n=%d\n' 'error n=%d\n' 'step mean %.2f%%\n' 'preStep SD %.2f%%\n' 'step z %.2f'], ...
    kernel.nCorrect, kernel.nError,  kernel.stepMean,  kernel.preStepSD, kernel.stepMeanZPreSD);
text(ax, 0.02, 0.98, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
  'FontSize', 8);

legend(ax, 'Location', 'southeast', 'FontSize', 7);

end
%% -------------------------------------------------------------------------
function plotSessionPerformance(ax, sessionRecords)

plot(ax, sessionRecords.meanCorrectStepNoise, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 4);

xlabel(ax, 'Session');
ylabel(ax, 'P(hit), noise step');
title(ax, 'Noise Step Performance by Session');
ylim(ax, [0.45 1.02]);
grid(ax, 'on');
box(ax, 'off');

end

%% -------------------------------------------------------------------------
function plotSessionKernelMeans(ax, sessionRecords)

hold(ax, 'on');

plot(ax, sessionRecords.kernelStepMean, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 4);
yline(ax, 0, 'k:');
xlabel(ax, 'Session');
ylabel(ax, 'Kernel Mean Over Step');
ytickformat(ax, 'percentage');
title(ax, 'Session Summed-Noise Kernel Means');
grid(ax, 'on');
box(ax, 'off');

end

function fig = plotIDQDirectionDiagnosticsSummary(acrossSummary)
% plotIDQDirectionDiagnosticsSummary
%
% Diagnostic page for direction-specific IDQ behavior and kernels.
% These diagnostics are not the primary readout analysis.

Dabs = acrossSummary.directionDiagnostics.absolute;
Daln = acrossSummary.directionDiagnostics.aligned;

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 11 8.5], ...
        'PaperOrientation', 'landscape', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title(tl, 'IDQ direction diagnostics', ...
  'Interpreter', 'none', ...
  'FontWeight', 'bold');

ax = nexttile(tl, 1);
plotDirectionBehavior(ax, Dabs.behavior);

ax = nexttile(tl, 2);
plotDirectionKernelSummary(ax, Dabs.kernels, ...
  'Absolute direction kernel step means');

ax = nexttile(tl, 3);
plotRectGainByDirection(ax, ...
    Dabs.rectGainByDirection);

ax = nexttile(tl, 4);
plotDirectionKernels(ax, Dabs.kernels, ...
  'Absolute changed-side kernels');

ax = nexttile(tl, 5);
plotDirectionKernels(ax, Daln.kernels, ...
  'Aligned changed-side kernels');

ax = nexttile(tl, 6);
plotDirectionDiagnosticsText(ax, acrossSummary);

end

%% -------------------------------------------------------------------------
function plotDirectionBehavior(ax, behavior)

bar(ax, behavior.pCorrect);
hold(ax, 'on');

yline(ax, 0.5, 'k:', 'HandleVisibility', 'off');

set(ax, ...
    'XTick', 1:height(behavior), ...
    'XTickLabel', behavior.directionLabel);

xtickangle(ax, 30);

ylim(ax, [0.45 1.0]);

ylabel(ax, 'P(correct)');
title(ax, 'Behavior by physical drift direction', ...
    'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');

for i = 1:height(behavior)
    text(ax, i, behavior.pCorrect(i), ...
        sprintf('n=%d', behavior.nTrials(i)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 7);
end

end

%% -------------------------------------------------------------------------
function plotDirectionKernelSummary(ax, kernels, plotTitle)

S = kernels.summaryTable;

bar(ax, S.stepMean);
hold(ax, 'on');

yline(ax, 0, 'k:', 'HandleVisibility', 'off');

set(ax, ...
    'XTick', 1:height(S), ...
    'XTickLabel', S.directionLabel);

xtickangle(ax, 30);

ylabel(ax, 'Step mean kernel');
title(ax, plotTitle, ...
    'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');

for i = 1:height(S)
    y = S.stepMean(i);

    if y >= 0
        va = 'bottom';
    else
        va = 'top';
    end

    text(ax, i, y, ...
        sprintf('%.3g', y), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', va, ...
        'FontSize', 7);
end
yl = ylim(ax);
m = max(abs(yl));
ylim(ax, [-m m]);
end

%% -------------------------------------------------------------------------
function plotDirectionKernels(ax, kernels, plotTitle)

hold(ax, 'on');

tMS = kernels.tMS;
stepFrames = kernels.stepFrames;

plot(ax, tMS, zeros(size(tMS)), 'k:', ...
    'HandleVisibility', 'off');

x1 = tMS(stepFrames(1));
x2 = tMS(stepFrames(end));

% Initial plot to set y limits.
for iDir = 1:numel(kernels.summaryTable.directionLabel)
    plot(ax, tMS, kernels.meanDiff(:, iDir), ...
        'LineWidth', 1.0, ...
        'DisplayName', char(kernels.summaryTable.directionLabel(iDir)));
end

yl = ylim(ax);

patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
    [0.9 0.9 0.9], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.35, ...
    'HandleVisibility', 'off');

% Replot on top of patch.
for iDir = 1:numel(kernels.summaryTable.directionLabel)
    plot(ax, tMS, kernels.meanDiff(:, iDir), ...
        'LineWidth', 1.0, ...
        'DisplayName', char(kernels.summaryTable.directionLabel(iDir)));
end

xlabel(ax, 'Time from trial start (ms)');
ylabel(ax, 'Correct - error noise');
title(ax, plotTitle, ...
    'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');

legend(ax, 'Location', 'best', 'FontSize', 7);

end

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
    'Direction diagnostics'
    ''
    'Behavior is grouped by physical drift direction.'
    'Aligned kernels are grouped by stream relative to drift.'
    ''
    sprintf('Best physical direction: %s, p=%.3f', ...
        B.directionLabel(bestDir), B.pCorrect(bestDir))
    sprintf('Worst physical direction: %s, p=%.3f', ...
        B.directionLabel(worstDir), B.pCorrect(worstDir))
    sprintf('Behavior range: %.3f', behaviorRange)
    ''
    'Absolute kernel step means:'
    sprintf('  %s: %.3g', Sabs.directionLabel(1), Sabs.stepMean(1))
    sprintf('  %s: %.3g', Sabs.directionLabel(2), Sabs.stepMean(2))
    sprintf('  %s: %.3g', Sabs.directionLabel(3), Sabs.stepMean(3))
    ''
    'Aligned kernel step means:'
    sprintf('  %s: %.3g', Saln.directionLabel(1), Saln.stepMean(1))
    sprintf('  %s: %.3g', Saln.directionLabel(2), Saln.stepMean(2))
    sprintf('  %s: %.3g', Saln.directionLabel(3), Saln.stepMean(3))
    ''
    'Diagnostic only; primary gain remains collapsed across directions.'
    };

text(ax, 0, 1, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', ...
    'FontSize', 8);

end

%% -------------------------------------------------------------------------
function plotRectGainByDirection(ax, rectGainByDirection)

S = rectGainByDirection.summaryTable;

hold(ax, 'on');

x = 1:height(S);

bar(ax, x, S.gain);

for i = 1:height(S)
  plot(ax, [x(i) x(i)], [S.CI95Low(i) S.CI95High(i)], ...
    'k-', 'LineWidth', 1.3);
end

plot(ax, x, S.gain, 'ko', ...
  'MarkerFaceColor', 'k', ...
  'MarkerSize', 5);

yline(ax, rectGainByDirection.flatPrediction, 'k--', ...
  'DisplayName', 'flat prediction');

yline(ax, 0, 'k:', ...
  'HandleVisibility', 'off');

set(ax, ...
  'XTick', x, ...
  'XTickLabel', S.directionLabel);

xtickangle(ax, 30);

ylabel(ax, 'Rect noise gain');
title(ax, 'Rectangular gain by drift direction', ...
  'Interpreter', 'none');

grid(ax, 'on');
box(ax, 'off');

yl = ylim(ax);
ylim(ax, [min([yl(1), min(S.CI95Low)-0.1, -0.1]), ...
  max([yl(2), rectGainByDirection.flatPrediction+0.1])]);

for i = 1:height(S)
  text(ax, x(i), S.CI95High(i), ...
    sprintf('zF %.1f', S.zVsFlat(i)), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 7);
end

end