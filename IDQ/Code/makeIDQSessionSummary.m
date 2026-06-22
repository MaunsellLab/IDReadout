function out = makeIDQSessionSummary(varargin)
% makeIDQSessionSummary
%
% Create daily single-session IDQ summary PDFs from processed session files.
%
% This function reads only:
%
%   Data/ProcessedSessions/ relative to: domainFolder(mfilename('fullpath'))
%
% Required variables in each processed session file:
%   header
%   sessionHeader
%   trialData
%   noiseBySideDir
%
% noiseBySideDir:
%   2 sides x 3 directions x nFrames x nTrials
%
% Conventions:
%   side dim 1 = RF, 2 = Opp
%   dir dim 1:3 matches sessionHeader.dirsDeg
%   trialData.sideIndex is changed side
%   trialData.chosenSideIndex is subject-chosen side
%   trialData.dirIndex is step direction index
%   trialData.correct is logical
%   trialData.stepCoh is step coherence
%   trialData.hasStepNoise is definitive inclusion flag for kernel analyses
%
% Example:
%   makeIDQSessionSummary
%   makeIDQSessionSummary('Replace', true)
%   makeIDQSessionSummary('ProcessedFile', 'IDQ_Animal_20260620.mat')
%
% Name-value inputs:
%   'ProcessedFile' : string/char/cellstr; default processes all .mat files
%   'Replace'       : logical; default false
%   'CloseFigures'  : logical; default true
%
% Output:
%   out is a struct array with processed file, PDF path, summary path, status.

p = inputParser;
p.addParameter('ProcessedFile', [], @(x) isempty(x) || ischar(x) || isstring(x) || iscellstr(x));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('CloseFigures', false, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

idqFolder = domainFolder(mfilename('fullpath'));
processedFolder = fullfile(idqFolder, 'Data', 'ProcessedSessions');
if ~isfolder(processedFolder)
    error('makeIDQSessionSummary:MissingProcessedFolder', ...
        'Processed session folder does not exist: %s', processedFolder);
end

dataOutFolder = validFolder(fullfile(idqFolder, 'Data', 'SessionSummaries'));
plotOutFolder = validFolder(fullfile(idqFolder, 'Plots', 'SessionSummaries'));
processedFiles = resolveProcessedFiles(processedFolder, opts.ProcessedFile);
out = repmat(struct('processedFile', '', 'pdfFile', '', 'summaryFile', '', 'status', '', ...
    'message', ''), 0, 1);

for iFile = 1:numel(processedFiles)
    processedFile = processedFiles{iFile};
    [~, baseName] = fileparts(processedFile);
    summaryFile = fullfile(dataOutFolder, [baseName '_sessionSummary.mat']);
    pdfFile = fullfile(plotOutFolder, [baseName '_sessionSummary.pdf']);

    thisOut = struct('processedFile', processedFile, ...
        'pdfFile', pdfFile, 'summaryFile', summaryFile, 'status', '', 'message', '');
    if ~opts.Replace && isOutputCurrent(processedFile, pdfFile, summaryFile)
        thisOut.status = 'skipped_current';
        thisOut.message = 'Existing PDF and summary MAT are newer than processed file.';
        out(end+1, 1) = thisOut; %#ok<AGROW>
        continue
    end

    load(processedFile, 'header', 'sessionHeader', 'trialData', 'noiseBySideDir');
    sessionAnalysis = loadOrMakeIDQSessionAnalysis(processedFile, sessionHeader);
    sessionSummary = computeSessionSummary(header, sessionHeader, sessionAnalysis, trialData, noiseBySideDir, ...
                      processedFile);
    fig = plotSessionSummaryFigure(sessionSummary);
    exportgraphics(fig, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');

    save(summaryFile, 'sessionSummary');

    if opts.CloseFigures
        close(fig);
    end

    thisOut.status = 'created';
    thisOut.message = 'Created session summary PDF and MAT.';
    out(end+1,1) = thisOut; %#ok<AGROW>
end

end

% -------------------------------------------------------------------------
function processedFiles = resolveProcessedFiles(processedFolder, processedFileOpt)

if isempty(processedFileOpt)
    d = dir(fullfile(processedFolder, '*.mat'));
    processedFiles = fullfile({d.folder}, {d.name});
    processedFiles = processedFiles(:);
    return
end

if ischar(processedFileOpt) || isstring(processedFileOpt)
    processedFileOpt = cellstr(processedFileOpt);
end
processedFiles = cell(size(processedFileOpt));
for i = 1:numel(processedFileOpt)
    f = char(processedFileOpt{i});
    candidate = fullfile(processedFolder, f);
    if isfile(candidate)
      processedFiles{i} = candidate;
    else
      fprintf('makeIDQSessionSummary:MissingProcessedFile: Cannot find processed file: %s', f);
    end
end
processedFiles = processedFiles(:);
end

% -------------------------------------------------------------------------
function tf = isOutputCurrent(processedFile, pdfFile, summaryFile)

if ~isfile(pdfFile) || ~isfile(summaryFile)
    tf = false;
    return
end

dProcessed = dir(processedFile);
dPdf = dir(pdfFile);
dSummary = dir(summaryFile);

tf = dPdf.datenum >= dProcessed.datenum && ...
     dSummary.datenum >= dProcessed.datenum;

end

% -------------------------------------------------------------------------
function sessionSummary = computeSessionSummary( ...
              header, sessionHeader, sessionAnalysis, trialData, noiseBySideDir, processedFile)

nTrials = size(noiseBySideDir, 4);
nFrames = size(noiseBySideDir, 3);

hasStepNoise = logical(trialData.hasStepNoise(:));
correct = logical(trialData.correct(:));
stepCoh = abs(double(trialData.stepCoh(:)));
dirIndex = double(trialData.dirIndex(:));
% sideIndex = double(trialData.sideIndex(:));
% stepSignIndex = double(trialData.stepSignIndex(:));

sessionSummary = struct;
sessionSummary.processedFile = processedFile;
sessionSummary.header = header;
sessionSummary.sessionHeader = sessionHeader;
sessionSummary.nTrials = nTrials;
sessionSummary.nFrames = nFrames;
sessionSummary.nValidTrials = numel(trialData.validIdx);
sessionSummary.nStepNoiseTrials = sum(hasStepNoise);
sessionSummary.nNoStepNoiseTrials = sum(~hasStepNoise);
sessionSummary.dirsDeg = double(sessionHeader.dirsDeg(:)');
sessionSummary.stepCohValues = unique(stepCoh(:))';

sessionSummary.psychometric.overall = summarizePsychometric(stepCoh, correct, hasStepNoise);
for iDir = 1:3
    idx = dirIndex == iDir;
    sessionSummary.psychometric.byDir(iDir) = summarizePsychometric(stepCoh(idx), correct(idx), hasStepNoise(idx)); 
end
sessionSummary.kernel = makeKernelSummaryFromSessionAnalysis(sessionAnalysis);
sessionSummary.plotLimits = choosePsychometricLimits(sessionSummary.psychometric);

%summarize biases
sessionSummary.bias = struct;
sessionSummary.bias.sideBias = struct;
sessionSummary.bias.sideBias.present = true;
sessionSummary.bias.sideBias.sourceField = 'sideBias';
sessionSummary.bias.sideBias.n = numel(trialData.sideBias);
sessionSummary.bias.sideBias.median = median(trialData.sideBias);
sessionSummary.bias.sideBias.iqr = prctile(trialData.sideBias, [25 75]);
sessionSummary.bias.sideBias.min = min(trialData.sideBias);
sessionSummary.bias.sideBias.max = max(trialData.sideBias);
sessionSummary.bias.sideBias.shape = [];

dirBiasMag = cellfun(@(s) s.magnitude, trialData.dirBias(:));
dirBiasMag(1:find(dirBiasMag ~= -1, 1, 'first')-1) = NaN;
sessionSummary.bias.dirBias = struct;
sessionSummary.bias.dirBias.present = true;
sessionSummary.bias.dirBias.sourceField = 'dirBias';
sessionSummary.bias.dirBias.n = numel(dirBiasMag);
sessionSummary.bias.dirBias.median = median(dirBiasMag, 'omitnan');
sessionSummary.bias.dirBias.iqr = prctile(dirBiasMag, [25 75]);
sessionSummary.bias.dirBias.min = min(dirBiasMag);
sessionSummary.bias.dirBias.max = max(dirBiasMag);
sessionSummary.bias.dirBias.shape = [];

end

% -------------------------------------------------------------------------
function psych = summarizePsychometric(stepCoh, correct, hasStepNoise)

stepCoh = stepCoh(:);
correct = logical(correct(:));
hasStepNoise = logical(hasStepNoise(:));

if isempty(stepCoh)
  emptyGroup = struct( ...
    'label', '', 'hasStepNoiseValue', [], 'coh', [],  'pCorrect', [],  'ciLow', [], 'ciHigh', []);
  psych = struct;
  psych.groups = repmat(emptyGroup, 1, 2);
  psych.combined = struct('coh', [], 'n', [], 'pCorrect', [], 'ciLow', [], 'ciHigh', []);
  psych.fit = struct('alpha', NaN, 'beta', NaN, 'xGrid', [], 'pFit', [], 'exitflag', NaN);
  return
end

cohVals = unique(stepCoh(:))';

groups = struct('label', {}, 'hasStepNoiseValue', {}, 'coh', {}, 'n', {}, 'pCorrect', {}, 'ciLow', {}, 'ciHigh', {});
groups(1) = summarizePsychometricGroup(stepCoh, correct, hasStepNoise, true, 'step noise', cohVals);
groups(2) = summarizePsychometricGroup(stepCoh, correct, hasStepNoise, false, 'no step noise', cohVals);
combined = summarizePsychometricCombined(stepCoh, correct, cohVals);

% fit = fitTwoAFCLogistic(stepCoh, correct);
fit = fitTwoAFCWeibull(stepCoh, correct);

psych = struct;
psych.groups = groups;
psych.combined = combined;
psych.fit = fit;

end

% -------------------------------------------------------------------------
function g = summarizePsychometricGroup(stepCoh, correct, hasStepNoise, hasNoiseValue, label, cohVals)

g = struct;
g.label = label;
g.hasStepNoiseValue = hasNoiseValue;
g.coh = cohVals;
g.n = zeros(size(cohVals));
g.pCorrect = nan(size(cohVals));
g.ciLow = nan(size(cohVals));
g.ciHigh = nan(size(cohVals));

for i = 1:numel(cohVals)
    idx = stepCoh == cohVals(i) & hasStepNoise == hasNoiseValue;
    g.n(i) = sum(idx);

    if g.n(i) > 0
        k = sum(correct(idx));
        g.pCorrect(i) = k / g.n(i);
        [g.ciLow(i), g.ciHigh(i)] = wilsonBinomialCI(k, g.n(i));
    end
end

end

% -------------------------------------------------------------------------
function c = summarizePsychometricCombined(stepCoh, correct, cohVals)

c = struct;
c.coh = cohVals;
c.n = zeros(size(cohVals));
c.pCorrect = nan(size(cohVals));
c.ciLow = nan(size(cohVals));
c.ciHigh = nan(size(cohVals));

for i = 1:numel(cohVals)
    idx = stepCoh == cohVals(i);
    c.n(i) = sum(idx);

    if c.n(i) > 0
        k = sum(correct(idx));
        c.pCorrect(i) = k / c.n(i);
        [c.ciLow(i), c.ciHigh(i)] = wilsonBinomialCI(k, c.n(i));
    end
end

end

%% -------------------------------------------------------------------------
function fit = fitTwoAFCWeibull(x, correct)

x = double(x(:));
y = double(correct(:));

valid = isfinite(x) & isfinite(y) & x >= 0;
x = x(valid);
y = y(valid);

fixedBetaShape = 3.0;
maxLapse = 0.05;

fit = struct( ...
  'alpha', NaN, ...
  'beta', fixedBetaShape, ...
  'lapse', NaN, ...
  'xGrid', [], ...
  'pFit', [], ...
  'exitflag', NaN);

if numel(unique(x)) < 2 || numel(unique(y)) < 2
  if ~isempty(x)
    fit.xGrid = linspace(0, max(x), 200);
    fit.pFit = nan(size(fit.xGrid));
  end
  return
end

xScale = max(abs(x));
if xScale == 0
  xScale = 1;
end
xs = x ./ xScale;

% Parameters are fit in unconstrained form:
%   rawAlpha -> alpha > 0
%   rawLapse -> 0 <= lapse <= maxLapse
%
% betaShape is fixed at 3.0 for stable session-level fits.

alpha0 = median(xs(xs > 0));
if isempty(alpha0) || ~isfinite(alpha0) || alpha0 <= 0
  alpha0 = 0.5;
end

lapse0 = 0.02;
lapse0 = min(max(lapse0, eps), maxLapse - eps);

rawAlpha0 = log(alpha0);
rawLapse0 = log(lapse0 / (maxLapse - lapse0));

b0 = [rawAlpha0 rawLapse0];

obj = @(b) negLogLikelihoodTwoAFCWeibullFixedBeta(b, xs, y, fixedBetaShape, maxLapse);

opts = optimset( ...
  'Display', 'off', ...
  'MaxFunEvals', 1e5, ...
  'MaxIter', 1e5, ...
  'TolX', 1e-8, ...
  'TolFun', 1e-8);

[bhat, ~, exitflag] = fminsearch(obj, b0, opts);

[alphaS, lapse] = unpackWeibullParamsFixedBeta(bhat, maxLapse);

fit.alpha = alphaS * xScale;
fit.beta = fixedBetaShape;
fit.lapse = lapse;
fit.exitflag = exitflag;

fit.xGrid = linspace(0, max(x), 200);
fit.pFit = twoAFCWeibull(fit.xGrid ./ xScale, alphaS, fixedBetaShape, lapse);

end

%% -------------------------------------------------------------------------
function nll = negLogLikelihoodTwoAFCWeibullFixedBeta(b, x, y, betaShape, maxLapse)

[alpha, lapse] = unpackWeibullParamsFixedBeta(b, maxLapse);

p = twoAFCWeibull(x, alpha, betaShape, lapse);

epsP = 1e-12;
p = min(max(p, epsP), 1 - epsP);

nll = -sum(y .* log(p) + (1 - y) .* log(1 - p));

end

%% -------------------------------------------------------------------------
function [alpha, lapse] = unpackWeibullParamsFixedBeta(b, maxLapse)

alpha = exp(b(1));              % alpha > 0
lapse = maxLapse ./ (1 + exp(-b(2)));

end

%% -------------------------------------------------------------------------
function p = twoAFCWeibull(x, alpha, betaShape, lapse)

gamma = 0.5;

x = max(x, 0);

p = gamma + (1 - gamma - lapse) .* ...
  (1 - exp(-(x ./ alpha) .^ betaShape));

end

% -------------------------------------------------------------------------
function [lo, hi] = wilsonBinomialCI(k, n)

if n == 0
    lo = NaN;
    hi = NaN;
    return
end

z = 1.96;
phat = k / n;
den = 1 + z^2 / n;
center = (phat + z^2 / (2*n)) / den;
halfWidth = z * sqrt((phat*(1 - phat) + z^2/(4*n)) / n) / den;

lo = max(0, center - halfWidth);
hi = min(1, center + halfWidth);

end

% % -------------------------------------------------------------------------
% function kernel = computeCombinedChangeSideKernel(noiseBySideDir, sideIndex, stepSignIndex, correct, hasStepNoise)
% 
% nTrials = size(noiseBySideDir, 4);
% nFrames = size(noiseBySideDir, 3);
% 
% if numel(sideIndex) ~= nTrials || ...
%    numel(stepSignIndex) ~= nTrials || ...
%    numel(correct) ~= nTrials || ...
%    numel(hasStepNoise) ~= nTrials
%     error('makeIDQSessionSummary:KernelInputSizeMismatch', ...
%         'Kernel input vectors must match number of trials.');
% end
% 
% idxUse = logical(hasStepNoise(:));
% if ~any(idxUse)
%     kernel = struct( ...
%         'meanDiff', nan(1, nFrames), ...
%         'correctMean', nan(1, nFrames), ...
%         'errorMean', nan(1, nFrames), ...
%         'nCorrect', 0, ...
%         'nError', 0, ...
%         'nTrials', 0, ...
%         'timeIndex', 1:nFrames);
%     return
% end
% 
% noiseChangedMeanDir = nan(sum(idxUse), nFrames);
% useTrialInds = find(idxUse);
% 
% for i = 1:numel(useTrialInds)
%     tr = useTrialInds(i);
%     side = sideIndex(tr);
%     raw = squeeze(mean(noiseBySideDir(side, :, :, tr), 2));
%     if iscolumn(raw)
%         raw = raw';
%     end
%     noiseChangedMeanDir(i, :) = double(raw);
% end
% 
% correctUse = correct(idxUse);
% 
% correctNoise = noiseChangedMeanDir(correctUse, :);
% errorNoise = noiseChangedMeanDir(~correctUse, :);
% 
% kernel = struct;
% kernel.correctMean = mean(correctNoise, 1, 'omitnan');
% kernel.errorMean = mean(errorNoise, 1, 'omitnan');
% kernel.meanDiff = kernel.correctMean - kernel.errorMean;
% kernel.nCorrect = size(correctNoise, 1);
% kernel.nError = size(errorNoise, 1);
% kernel.nTrials = size(noiseChangedMeanDir, 1);
% kernel.timeIndex = 1:nFrames;
% 
% end

% -------------------------------------------------------------------------
function plotLimits = choosePsychometricLimits(psychometric)

allCoh = [];
allCiLow = [];

names = fieldnames(psychometric);
for iName = 1:numel(names)
    value = psychometric.(names{iName});

    if strcmp(names{iName}, 'overall')
        allCoh = [allCoh value.combined.coh]; %#ok<AGROW>
        allCiLow = [allCiLow value.combined.ciLow]; %#ok<AGROW>
    elseif strcmp(names{iName}, 'byDir')
        for iDir = 1:numel(value)
            allCoh = [allCoh value(iDir).combined.coh]; %#ok<AGROW>
            allCiLow = [allCiLow value(iDir).combined.ciLow]; %#ok<AGROW>
        end
    end
end

allCoh = allCoh(isfinite(allCoh));
allCiLow = allCiLow(isfinite(allCiLow));

if isempty(allCoh)
    xMax = 1;
else
    xMax = max(allCoh);
    if xMax == 0
        xMax = 1;
    end
end

if isempty(allCiLow)
    yMin = 0.45;
else
    yMin = min(allCiLow) - 0.05;
    yMin = max(0, floor(yMin * 20) / 20);
end

plotLimits = struct;
plotLimits.xLim = [0, xMax * 1.05];
plotLimits.yLim = [yMin, 1];

end

% -------------------------------------------------------------------------
function fig = plotSessionSummaryFigure(sessionSummary)

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 14 8]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
[~, baseName] = fileparts(sessionSummary.processedFile);
ttl = sprintf('IDQ session summary: %s', baseName);
title(tl, ttl, 'Interpreter', 'none', 'FontWeight', 'bold');

axText = nexttile(tl, 1);
plotTextSummary(axText, sessionSummary);

axPsych = nexttile(tl, 2);
plotPsychometricPanel(axPsych, sessionSummary.psychometric.overall, ...
    'Overall', sessionSummary.plotLimits);

axKernel = nexttile(tl, 3);
plotCombinedKernelPanel(axKernel, sessionSummary.kernel);

for iDir = 1:3
    ax = nexttile(tl, 3 + iDir);
    label = sprintf('%.0f°', sessionSummary.dirsDeg(iDir));
    plotPsychometricPanel(ax, sessionSummary.psychometric.byDir(iDir), label, sessionSummary.plotLimits);
end

end

% -------------------------------------------------------------------------
function plotTextSummary(ax, sessionSummary)

axis(ax, 'off');

lines = {};

[~, baseName, ~] = fileparts(sessionSummary.processedFile);
lines{end+1} = sprintf('%s', baseName);
lines{end+1} = sprintf('Trials: %d', sessionSummary.nTrials);
lines{end+1} = sprintf('Trials with Step Noise: %d', sessionSummary.nStepNoiseTrials);
lines{end+1} = sprintf('Directions: %d°, %d°, %d°', sessionSummary.dirsDeg(1:3));
cohStr = 'Step coherences:';
for c = 1:numel(sessionSummary.stepCohValues)
  cohStr = sprintf('%s %ld%%', cohStr, sessionSummary.stepCohValues(c));
end
lines{end+1} = cohStr;
lines{end+1} = '';
lines{end+1} = 'Biases';

sb = sessionSummary.bias.sideBias;
if sb.present && sb.n > 0
    lines{end+1} = sprintf('  Side Bias: median %.3g, IQR [%.3g %.3g]', sb.median, sb.iqr(1), sb.iqr(2));
else
    lines{end+1} = 'sideBias: not found';
end

db = sessionSummary.bias.dirBias;
if db.present && db.n > 0
    lines{end+1} = sprintf('  Direction bias: median %.3g, IQR [%.3g %.3g]', db.median, db.iqr(1), db.iqr(2));
else
    lines{end+1} = 'direction bias: not found';
end

lines{end+1} = '';
lines{end+1} = 'Kernel Evidence: Drift Noise - Non-Drift Noise';
text(ax, 0, 1, strjoin(lines, newline), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',  'FontName', 'Helvetica', ...
    'FontSize', 9, 'Interpreter', 'none');

end

% -------------------------------------------------------------------------
function plotPsychometricPanel(ax, psych, panelTitle, plotLimits)

hold(ax, 'on');

plot(ax, [plotLimits.xLim(1) plotLimits.xLim(2)], [0.5 0.5], ':', ...
    'HandleVisibility', 'off');

markerList = {'o', 's'};

for iGroup = 1:numel(psych.groups)
    g = psych.groups(iGroup);
    idx = g.n > 0 & isfinite(g.pCorrect);

    if ~any(idx)
        continue
    end

    x = g.coh(idx);
    y = g.pCorrect(idx);
    lo = g.ciLow(idx);
    hi = g.ciHigh(idx);

    errLow = y - lo;
    errHigh = hi - y;

    errorbar(ax, x, y, errLow, errHigh, markerList{iGroup}, ...
        'LineStyle', 'none', ...
        'MarkerSize', 5, ...
        'DisplayName', g.label);
    n = g.n(idx);
    for i = 1:numel(x)
      if iGroup == 1
        dy = 0.025;
      else
        dy = -0.035;
      end
      text(ax, x(i), y(i) + dy, sprintf('%d', n(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 7, 'Clipping', 'on');
    end
end
if isfield(psych, 'fit') && ~isempty(psych.fit.xGrid) && any(isfinite(psych.fit.pFit))
    plot(ax, psych.fit.xGrid, psych.fit.pFit, '-', ...
        'LineWidth', 1.2, ...
        'DisplayName', 'combined fit');
end
addWeibullFitTextBox(gca, psych.fit);

title(ax, panelTitle, 'Interpreter', 'none');
xlabel(ax, 'Coherence Step (from Pre-Step, %)');
ylabel(ax, 'Proportion Correct');
xlim(ax, plotLimits.xLim);
ylim(ax, plotLimits.yLim);
grid(ax, 'on');
box(ax, 'off');

legend(ax, 'Location', 'southwest', 'FontSize', 7);
end

%% -------------------------------------------------------------------------
function addWeibullFitTextBox(ax, fit)

if nargin < 1 || isempty(ax) || ~isgraphics(ax, 'axes')
  ax = gca;
end

if nargin < 2 || isempty(fit) || ~isstruct(fit)
  return
end

requiredFields = {'alpha', 'beta', 'lapse'};
if ~all(isfield(fit, requiredFields))
  return
end

if ~isfinite(fit.alpha) || ~isfinite(fit.beta) || ~isfinite(fit.lapse)
  return
end

txt = sprintf([ ...
  '2AFC Weibull fit\n', ...
  'P = 0.5 + (0.5 - \\lambda)(1 - e^{-(x/\\alpha)^\\beta})\n', ...
  '\\alpha = %.3g %% coh\n', ...
  '\\beta = %.3g\n', ...
  '\\lambda = %.3g'], ...
  fit.alpha, fit.beta, fit.lapse);

text(ax, 0.97, 0.03, txt, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
  'FontSize', 8,  'Interpreter', 'tex', 'BackgroundColor', 'w', 'EdgeColor', [0.4 0.4 0.4],  'Margin', 2);

end

%% -------------------------------------------------------------------------
function plotCombinedKernelPanel(ax, kernel)

hold(ax, 'on');

plot(ax, kernel.timeIndex, zeros(size(kernel.timeIndex)), ':', ...
    'HandleVisibility', 'off');

plot(ax, kernel.timeIndex, kernel.meanDiff, '-', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'correct - error');

title(ax, sprintf('Changed-side signed-noise kernel, n=%d', kernel.nTrials), ...
    'Interpreter', 'none');
xlabel(ax, 'Frame');
ylabel(ax, 'Signed noise, correct - error');
grid(ax, 'on');
box(ax, 'off');

legend(ax, 'Location', 'best', 'FontSize', 7);

txt = sprintf('hit n=%d\nmiss n=%d', kernel.nCorrect, kernel.nError);
text(ax, 0.02, 0.98, txt, ...
    'Units', 'normalized', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'FontSize', 8);

end

%% ------------------------------------------------------------------------
function sessionAnalysis = loadOrMakeIDQSessionAnalysis(processedFile, sessionHeader)

domainPath = domainFolder(mfilename('fullpath'));
analysisFolder = fullfile(domainPath, 'Data', 'SessionAnalysis');
sessionName = sessionHeader.fileName;
sessionName = erase(sessionName, '.mat');
sessionName = erase(sessionName, '.dat');
analysisFile = fullfile(analysisFolder, sprintf('%s_sessionAnalysis.mat', sessionName));
if ~isfile(analysisFile)
  makeIDQSessionAnalysis(processedFile);
end
load(analysisFile, 'sessionAnalysis');
end

%% -------------------------------------------------------------------------
function kernel = makeKernelSummaryFromSessionAnalysis(sessionAnalysis)

kernel = struct();

kernel.correctMean = sessionAnalysis.signedNoiseKernel.meanCorrect';
kernel.errorMean = sessionAnalysis.signedNoiseKernel.meanError';
kernel.meanDiff = sessionAnalysis.signedNoiseKernel.kernel';

kernel.nCorrect = sessionAnalysis.signedNoiseKernel.nCorrect;
kernel.nError = sessionAnalysis.signedNoiseKernel.nError;
kernel.nTrials = sessionAnalysis.signedNoiseKernel.nTrials;

kernel.timeIndex = 1:numel(sessionAnalysis.tMS);

end