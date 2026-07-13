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
% fprintf('plotIDQSummary: loading %d sessionAnalysis files\n', numel(files));

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

targetPerformance = 0.75;
noiseGain = fitIDQNoiseGain(trialTable, psych.sessionFits, psych.alignedWeibull, targetPerformance);
directionDiagnostics = computeIDQDirectionDiagnostics(sessionAnalyses, trialTable);

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

acrossSummary.noiseGain = noiseGain;
acrossSummary.directionDiagnostics = directionDiagnostics;

summaryMatFile = fullfile(summaryFolder, 'IDQ_AcrossSessionSummary.mat');
save(summaryMatFile, 'acrossSummary', '-v7.3');
% fprintf('  saved %s\n', summaryMatFile);

fig = plotSummary(acrossSummary);
exportgraphics(fig, fullfile(plotFolder, 'IDQSessionSummary.pdf'), 'ContentType', 'vector');

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

fit = psych.alignedWeibull;

% fprintf('lapse       = %.8f\n', fit.lapse);
% fprintf('lapseRaw    = %.8f\n', fit.thetaHat(3));
% fprintf('threshold   = %.6f\n', fit.threshold);
% fprintf('beta        = %.6f\n', fit.betaWeibull);
% fprintf('NLL         = %.6f\n', fit.negLogLikelihood);
% fprintf('exitflag    = %d\n', fit.exitflag);
% disp(fit.hessian)

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
  k = SA.sumNoiseKernel.kernel;
  kStep = k(SA.stepFrames);
  kernelStepIntegral(iSession) = sum(kStep, 'omitnan');
  kernelStepMean(iSession) = mean(kStep, 'omitnan');
  kernelStepPeak(iSession) = max(abs(kStep), [], 'omitnan');
end

sessionRecords = table(fileName, nTrials, nStepNoiseTrials, noisyStepCoh, meanCorrect, ...
    meanCorrectStepNoise, meanRectNoise, stdRectNoise, kernelStepIntegral, kernelStepMean, kernelStepPeak);

end

%% -------------------------------------------------------------------------
function fig = plotSummary(acrossSummary)

fig = figure(500);
set(fig, 'Color', 'w', 'WindowStyle', 'docked');
tl = tiledlayout(fig, 4, 4, 'TileSpacing', 'loose', 'Padding', 'compact');
title(tl, sprintf('IDQ Summary (%d sessions)', acrossSummary.nSessions), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

textAx = nexttile(tl, 12);
plotTextSummary(textAx, acrossSummary);
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
  stepPatch(d) = plotKernel(kernelAx(d), acrossSummary.absKernels{d}, absColor);
end

kernelAx(4) = nexttile(tl, 8);
stepPatch(4) = plotKernel(kernelAx(4), acrossSummary.kernel, allColor);

for d = 1:3
  kernelAx(d + 4) = nexttile(tl, d+  8);
  stepPatch(d + 4) = plotKernel(kernelAx(d + 4), acrossSummary.relKernels{d}, relColor, d == 1, false);
end

yl = cell2mat(get(kernelAx, 'YLim'));
sharedYLim = [min(yl(:,1)), max(yl(:,2))];
set(kernelAx, 'YLim', sharedYLim);
for k = 1:7
  stepPatch(k).YData = sharedYLim([1 1 2 2]);
end

% Use a common y-axis scale for all gain plots.
gainYLim = getNoiseGainYLim(acrossSummary.noiseGain);
plotGainBars(nexttile(tl, 13), acrossSummary.noiseGain.driftRelative, 'Drift-Rel. Gains', gainYLim, relColor);
plotGainBars(nexttile(tl, 14), acrossSummary.noiseGain.absolute, 'Abs. Direction Gains', gainYLim, absColor);
plotGainBars(nexttile(tl, 15), acrossSummary.noiseGain.combined, 'Combined Gain', gainYLim, allColor);
plotGainBars(nexttile(tl, 16), acrossSummary.noiseGain.driftNonDrift, 'Drift v. Non-Drift Gains', gainYLim, relColor);

end

%% -------------------------------------------------------------------------
function plotTextSummary(ax, acrossSummary)

axis(ax, 'off');
txt = {
    sprintf('Created: %s', string(acrossSummary.createdAt, 'dd/MM/yyyy HH:mm'))
    };

text(ax, 0, 1, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontName', 'Menlo', 'FontSize', 8);
end

%% -------------------------------------------------------------------------
function stepPatch = plotKernel(ax, kernel, plotColor, showLegend, showXLabel)

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
ylabel(ax, 'Change Side Kernel');
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