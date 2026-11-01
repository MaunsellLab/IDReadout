function results = makeIDRProbeOffsetThresholdFirstPass(fitPath, varargin)
%MAKEIDRPROBEOFFSETTHRESHOLDFIRSTPASS Descriptive IDR threshold-offset test.
%
% results = makeIDRProbeOffsetThresholdFirstPass(fitPath)
%
% This first-pass analysis deliberately avoids a mechanistic noise model. For
% each one-coherence probeSession it converts observed P(correct) to the 75%
% threshold implied by a fixed Weibull beta and lapse, then divides that value
% by the parent session's no-noise threshold. It also reports the simpler
% observed-minus-expected P(correct) residual.
%
% The default beta and lapse are the shared values stored in fitData.probeTable,
% keeping the conversion consistent with alphaParent. Parent sessions are the
% resampling unit. An offset-label permutation test uses only within-parent
% information (singleton parent sessions do not contribute to that test).
%
% Name-value options
%   'OutputFolder'      folder for PDF, MAT, and CSV outputs (default: fit folder)
%   'OutputStem'        output filename stem (default derived from input)
%   'Beta'              fixed Weibull beta (default: stored shared beta)
%   'Lapse'             fixed lapse (default: stored shared lapse)
%   'Combine179And180'  pool 180 deg with 179 deg (default true)
%   'NumBootstrap'      parent-session bootstrap samples (default 5000)
%   'NumPermutations'   within-parent permutations (default 10000)
%   'RandomSeed'        random seed (default 1)
%   'MakePlot'          create and export summary PDF (default true)
%   'SaveOutputs'       save MAT and CSV outputs (default true)

p = inputParser;
p.addRequired('fitPath', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('OutputFolder', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('OutputStem', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('Beta', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('Lapse', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0 && x < 0.25));
p.addParameter('Combine179And180', true, @(x) islogical(x) && isscalar(x));
p.addParameter('NumBootstrap', 5000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1) == 0);
p.addParameter('NumPermutations', 10000, @(x) isnumeric(x) && isscalar(x) && x >= 0 && mod(x,1) == 0);
p.addParameter('RandomSeed', 1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
p.addParameter('MakePlot', true, @(x) islogical(x) && isscalar(x));
p.addParameter('SaveOutputs', true, @(x) islogical(x) && isscalar(x));
p.parse(fitPath, varargin{:});
opts = p.Results;

fitPath = char(fitPath);
if ~isfile(fitPath)
  error('makeIDRProbeOffsetThresholdFirstPass:MissingInput', ...
    'Input MAT file was not found: %s', fitPath);
end

S = load(fitPath, 'fitData');
if ~isfield(S, 'fitData') || ~isfield(S.fitData, 'probeTable')
  error('makeIDRProbeOffsetThresholdFirstPass:MissingProbeTable', ...
    'The input must contain fitData.probeTable.');
end
T = S.fitData.probeTable;
required = {'probeSessionName','parentSessionName','probeDirDeg','alphaParent', ...
  'betaShared','lapseShared','nIncTrials','nCoherenceLevels','minCoherence', ...
  'maxCoherence','fractionCorrect','fitIncluded'};
missing = required(~ismember(required, T.Properties.VariableNames));
if ~isempty(missing)
  error('makeIDRProbeOffsetThresholdFirstPass:MissingVariables', ...
    'probeTable is missing: %s', strjoin(missing, ', '));
end

valid = logical(T.fitIncluded) & isfinite(T.probeDirDeg) & ...
  isfinite(T.alphaParent) & T.alphaParent > 0 & ...
  isfinite(T.nIncTrials) & T.nIncTrials > 0 & ...
  isfinite(T.fractionCorrect) & T.nCoherenceLevels == 1 & ...
  isfinite(T.minCoherence) & isfinite(T.maxCoherence);
if any(logical(T.fitIncluded) & T.nCoherenceLevels ~= 1)
  warning('makeIDRProbeOffsetThresholdFirstPass:MultipleCoherences', ...
    'Excluding fitted probeSessions that do not contain exactly one coherence.');
end
T = T(valid,:);
if isempty(T)
  error('makeIDRProbeOffsetThresholdFirstPass:NoValidRows', ...
    'No valid one-coherence probeSessions were found.');
end

if opts.Combine179And180
  T.probeDirDeg(T.probeDirDeg == 180) = 179;
end

storedBeta = median(double(T.betaShared), 'omitnan');
storedLapse = median(double(T.lapseShared), 'omitnan');
if isempty(opts.Beta)
  beta = storedBeta;
else
  beta = double(opts.Beta);
end
if isempty(opts.Lapse)
  lapse = storedLapse;
else
  lapse = double(opts.Lapse);
end
if ~isfinite(beta) || beta <= 0 || ~isfinite(lapse) || lapse < 0 || lapse >= 0.25
  error('makeIDRProbeOffsetThresholdFirstPass:InvalidShape', ...
    'A finite beta > 0 and lapse in [0, 0.25) are required.');
end

betaValues = double(T.betaShared);
lapseValues = double(T.lapseShared);
betaRange = max(betaValues, [], 'omitnan') - min(betaValues, [], 'omitnan');
lapseRange = max(lapseValues, [], 'omitnan') - min(lapseValues, [], 'omitnan');
if betaRange > 1e-8 || lapseRange > 1e-8
  warning('makeIDRProbeOffsetThresholdFirstPass:NonconstantShape', ...
    'Stored beta/lapse are not constant; their medians are being used.');
end

coherence = (double(T.minCoherence) + double(T.maxCoherence)) ./ 2;
nTrials = double(T.nIncTrials);
nCorrect = round(double(T.fractionCorrect) .* nTrials);
pObserved = nCorrect ./ nTrials;

% Keep the inverse Weibull finite for pathological proportions. This normally
% has no effect on IDR data, whose observed proportions lie well inside bounds.
lower = 0.5;
upper = 1 - lapse;
clipMargin = 0.5 ./ (nTrials + 1);
pForInverse = min(max(pObserved, lower + clipMargin), upper - clipMargin);
zObserved = -log(1 - (pForInverse - 0.5) ./ (0.5 - lapse));
z75 = -log(1 - 0.25 ./ (0.5 - lapse));

equivalentC75 = coherence .* (z75 ./ zObserved) .^ (1 ./ beta);
parentC75 = double(T.alphaParent) .* z75 .^ (1 ./ beta);
thresholdRatio = equivalentC75 ./ parentC75;
logThresholdRatio = log(thresholdRatio);
pExpected = weibullPcorrectLocal(coherence, double(T.alphaParent), beta, lapse);
performanceResidual = pObserved - pExpected;

probeSessionTable = table(string(T.probeSessionName), string(T.parentSessionName), ...
  double(T.probeDirDeg), nTrials, nCorrect, coherence, pObserved, pExpected, ...
  performanceResidual, double(T.alphaParent), parentC75, equivalentC75, ...
  thresholdRatio, logThresholdRatio, ...
  'VariableNames', {'probeSessionName','parentSessionName','probeDirDeg', ...
  'nTrials','nCorrect','coherence','pObserved','pExpected','performanceResidual', ...
  'alphaParent','parentC75','equivalentC75','thresholdRatio','logThresholdRatio'});

oldRng = rng;
cleanupRng = onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(opts.RandomSeed, 'twister');

offsets = unique(probeSessionTable.probeDirDeg(:))';
nOffsets = numel(offsets);
parents = unique(probeSessionTable.parentSessionName, 'stable');
nParents = numel(parents);

bootThreshold = nan(opts.NumBootstrap, nOffsets);
bootResidual = nan(opts.NumBootstrap, nOffsets);
for iBoot = 1:opts.NumBootstrap
  sampledParents = parents(randi(nParents, nParents, 1));
  rowSets = cell(nParents,1);
  for iParent = 1:nParents
    rowSets{iParent} = find(probeSessionTable.parentSessionName == sampledParents(iParent));
  end
  sampledRows = vertcat(rowSets{:});
  for iOffset = 1:nOffsets
    m = probeSessionTable.probeDirDeg(sampledRows) == offsets(iOffset);
    if any(m)
      rows = sampledRows(m);
      bootThreshold(iBoot,iOffset) = median(probeSessionTable.thresholdRatio(rows), 'omitnan');
      bootResidual(iBoot,iOffset) = median(probeSessionTable.performanceResidual(rows), 'omitnan');
    end
  end
end

offsetSummary = table('Size', [nOffsets 10], ...
  'VariableTypes', repmat({'double'}, 1, 10), ...
  'VariableNames', {'probeDirDeg','nProbeSessions','nParents','medianThresholdRatio', ...
  'thresholdRatioCILow','thresholdRatioCIHigh','medianPerformanceResidual', ...
  'performanceResidualCILow','performanceResidualCIHigh','medianPObserved'});
for iOffset = 1:nOffsets
  m = probeSessionTable.probeDirDeg == offsets(iOffset);
  offsetSummary.probeDirDeg(iOffset) = offsets(iOffset);
  offsetSummary.nProbeSessions(iOffset) = sum(m);
  offsetSummary.nParents(iOffset) = numel(unique(probeSessionTable.parentSessionName(m)));
  offsetSummary.medianThresholdRatio(iOffset) = median(probeSessionTable.thresholdRatio(m), 'omitnan');
  offsetSummary.medianPerformanceResidual(iOffset) = median(probeSessionTable.performanceResidual(m), 'omitnan');
  offsetSummary.medianPObserved(iOffset) = median(probeSessionTable.pObserved(m), 'omitnan');
  if opts.NumBootstrap > 0
    ci = percentileLocal(bootThreshold(:,iOffset), [2.5 97.5]);
    offsetSummary.thresholdRatioCILow(iOffset) = ci(1);
    offsetSummary.thresholdRatioCIHigh(iOffset) = ci(2);
    ci = percentileLocal(bootResidual(:,iOffset), [2.5 97.5]);
    offsetSummary.performanceResidualCILow(iOffset) = ci(1);
    offsetSummary.performanceResidualCIHigh(iOffset) = ci(2);
  else
    offsetSummary{iOffset, {'thresholdRatioCILow','thresholdRatioCIHigh', ...
      'performanceResidualCILow','performanceResidualCIHigh'}} = nan(1,4);
  end
end

[omnibusStatistic, omnibusP, nInterleavedParents] = ...
  withinParentPermutationTest(probeSessionTable, opts.NumPermutations);

results = struct();
results.version = 1;
results.method = ['One-point equivalent 75% thresholds with fixed shared beta/lapse; ' ...
  'thresholds normalized by parent no-noise c75; parent-session bootstrap; ' ...
  'within-parent offset-label permutation test.'];
results.sourceFitPath = fitPath;
results.options = opts;
results.beta = beta;
results.lapse = lapse;
results.storedBeta = storedBeta;
results.storedLapse = storedLapse;
results.probeSessionTable = probeSessionTable;
results.offsetSummary = offsetSummary;
results.bootstrapThresholdRatio = bootThreshold;
results.bootstrapPerformanceResidual = bootResidual;
results.omnibusStatistic = omnibusStatistic;
results.omnibusPermutationP = omnibusP;
results.nParentSessions = nParents;
results.nInterleavedParentSessions = nInterleavedParents;
results.nProbeSessions = height(probeSessionTable);

[inputFolder, inputStem] = fileparts(fitPath);
if isempty(opts.OutputFolder)
  outputFolder = inputFolder;
else
  outputFolder = char(opts.OutputFolder);
end
if isempty(outputFolder)
  outputFolder = pwd;
end
if isempty(opts.OutputStem)
  outputStem = sprintf('%s_ThresholdFirstPass', inputStem);
else
  outputStem = char(opts.OutputStem);
end
if (opts.MakePlot || opts.SaveOutputs) && ~isfolder(outputFolder)
  mkdir(outputFolder);
end

results.outputPaths = struct('pdf', '', 'mat', '', 'probeSessionCSV', '', 'offsetSummaryCSV', '');
if opts.MakePlot
  pdfPath = fullfile(outputFolder, [outputStem '.pdf']);
  makeSummaryPlot(probeSessionTable, offsetSummary, beta, lapse, ...
    omnibusP, nInterleavedParents, opts, pdfPath);
  results.outputPaths.pdf = pdfPath;
end
if opts.SaveOutputs
  matPath = fullfile(outputFolder, [outputStem '.mat']);
  probeCSVPath = fullfile(outputFolder, [outputStem '_ProbeSessions.csv']);
  summaryCSVPath = fullfile(outputFolder, [outputStem '_OffsetSummary.csv']);
  save(matPath, 'results', '-v7.3');
  writetable(probeSessionTable, probeCSVPath);
  writetable(offsetSummary, summaryCSVPath);
  results.outputPaths.mat = matPath;
  results.outputPaths.probeSessionCSV = probeCSVPath;
  results.outputPaths.offsetSummaryCSV = summaryCSVPath;
  % Re-save so the MAT copy contains its own finalized output paths.
  save(matPath, 'results', '-v7.3');
end
end

% -------------------------------------------------------------------------
function p = weibullPcorrectLocal(coherence, alpha, beta, lapse)
p = 0.5 + (0.5 - lapse) .* (1 - exp(-((coherence ./ alpha) .^ beta)));
end

% -------------------------------------------------------------------------
function q = percentileLocal(x, pct)
x = sort(x(isfinite(x)));
if isempty(x)
  q = nan(size(pct));
  return;
end
if numel(x) == 1
  q = repmat(x, size(pct));
  return;
end
pos = 1 + (numel(x) - 1) .* pct ./ 100;
lo = floor(pos);
hi = ceil(pos);
w = pos - lo;
q = x(lo) .* (1 - w) + x(hi) .* w;
end

% -------------------------------------------------------------------------
function [statistic, pValue, nInterleaved] = withinParentPermutationTest(T, nPermutations)
% Reduction in RSS when offset terms are added to parent fixed effects.
y = T.logThresholdRatio;
[~,~,parentIndex] = unique(T.parentSessionName, 'stable');
[~,~,offsetIndex] = unique(T.probeDirDeg, 'sorted');
nParents = max(parentIndex);
nOffsets = max(offsetIndex);
nUniqueOffsets = accumarray(parentIndex, offsetIndex, [nParents 1], ...
  @(x) numel(unique(x)));
nInterleaved = sum(nUniqueOffsets > 1);

Xparent = full(sparse((1:height(T))', parentIndex, 1, height(T), nParents));
if nOffsets > 1
  Xoffset = full(sparse((1:height(T))', offsetIndex, 1, height(T), nOffsets));
  Xfull = [Xparent Xoffset(:,2:end)];
else
  Xfull = Xparent;
end
projectionDifference = Xfull * pinv(Xfull) - Xparent * pinv(Xparent);
projectionDifference = (projectionDifference + projectionDifference') ./ 2;
statistic = max(y' * projectionDifference * y, 0);

if nPermutations == 0 || nOffsets < 2 || nInterleaved == 0
  pValue = NaN;
  return;
end

permStatistics = zeros(nPermutations,1);
parentRows = accumarray(parentIndex, (1:height(T))', [nParents 1], @(x) {x});
for iPerm = 1:nPermutations
  yPerm = y;
  for iParent = 1:nParents
    rows = parentRows{iParent};
    if nUniqueOffsets(iParent) > 1
      yPerm(rows) = y(rows(randperm(numel(rows))));
    end
  end
  permStatistics(iPerm) = max(yPerm' * projectionDifference * yPerm, 0);
end
pValue = (1 + sum(permStatistics >= statistic - 1e-12)) ./ (nPermutations + 1);
end

% -------------------------------------------------------------------------
function makeSummaryPlot(T, S, beta, lapse, omnibusP, nInterleaved, opts, pdfPath)
offsets = S.probeDirDeg(:)';
nOffsets = numel(offsets);
colors = lines(max(nOffsets, 7));

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 12 8]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'IDR probe-offset threshold first pass', 'FontWeight', 'bold');

ax1 = nexttile(layout, 1);
hold(ax1, 'on');
for iOffset = 1:nOffsets
  m = T.probeDirDeg == offsets(iOffset);
  jitter = deterministicJitter(sum(m), 0.18);
  scatter(ax1, iOffset + jitter, T.thresholdRatio(m), 17, ...
    'MarkerFaceColor', colors(iOffset,:), 'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.28);
end
errorbar(ax1, 1:nOffsets, S.medianThresholdRatio, ...
  S.medianThresholdRatio - S.thresholdRatioCILow, ...
  S.thresholdRatioCIHigh - S.medianThresholdRatio, ...
  'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.4, 'CapSize', 7);
yline(ax1, 1, ':', 'parent no-noise threshold');
formatOffsetAxis(ax1, offsets);
ylabel(ax1, 'Equivalent c_{75} / parent c_{75}');
title(ax1, 'Normalized equivalent thresholds');
grid(ax1, 'on'); box(ax1, 'off');

ax2 = nexttile(layout, 2);
hold(ax2, 'on');
for iOffset = 1:nOffsets
  m = T.probeDirDeg == offsets(iOffset);
  jitter = deterministicJitter(sum(m), 0.18);
  scatter(ax2, iOffset + jitter, T.performanceResidual(m), 17, ...
    'MarkerFaceColor', colors(iOffset,:), 'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.28);
end
errorbar(ax2, 1:nOffsets, S.medianPerformanceResidual, ...
  S.medianPerformanceResidual - S.performanceResidualCILow, ...
  S.performanceResidualCIHigh - S.medianPerformanceResidual, ...
  'ko', 'MarkerFaceColor', 'k', 'LineWidth', 1.4, 'CapSize', 7);
yline(ax2, 0, ':');
formatOffsetAxis(ax2, offsets);
ylabel(ax2, 'Observed - expected P(correct)');
title(ax2, 'Performance residuals');
grid(ax2, 'on'); box(ax2, 'off');

ax3 = nexttile(layout, 3);
hold(ax3, 'on');
parents = unique(T.parentSessionName, 'stable');
for iParent = 1:numel(parents)
  m = T.parentSessionName == parents(iParent);
  [parentOffsets, order] = sort(T.probeDirDeg(m));
  if numel(unique(parentOffsets)) < 2
    continue;
  end
  [~, x] = ismember(parentOffsets, offsets);
  values = T.thresholdRatio(m);
  plot(ax3, x, values(order), '-', 'Color', [0.72 0.72 0.72], 'LineWidth', 0.6);
end
plot(ax3, 1:nOffsets, S.medianThresholdRatio, '-ko', ...
  'MarkerFaceColor', 'k', 'LineWidth', 1.7);
yline(ax3, 1, ':');
formatOffsetAxis(ax3, offsets);
ylabel(ax3, 'Equivalent c_{75} / parent c_{75}');
title(ax3, sprintf('Within-parent connections (n = %d parents)', nInterleaved));
grid(ax3, 'on'); box(ax3, 'off');

ax4 = nexttile(layout, 4);
axis(ax4, 'off');
if isnan(omnibusP)
  pText = 'not computed';
else
  pText = sprintf('%.4g', omnibusP);
end
linesText = {
  sprintf('ProbeSessions: %d', height(T))
  sprintf('Parent sessions: %d', numel(unique(T.parentSessionName)))
  sprintf('Interleaved parent sessions: %d', nInterleaved)
  sprintf('Fixed beta: %.4f', beta)
  sprintf('Fixed lapse: %.5f', lapse)
  sprintf('Parent bootstrap samples: %d', opts.NumBootstrap)
  sprintf('Within-parent permutation p: %s', pText)
  ''
  'Offset     probeSessions     parents'
  };
for iOffset = 1:height(S)
  linesText{end+1,1} = sprintf('%5g%15d%12d', S.probeDirDeg(iOffset), ...
    S.nProbeSessions(iOffset), S.nParents(iOffset)); %#ok<AGROW>
end
text(ax4, 0.02, 0.98, strjoin(linesText, newline), ...
  'Units', 'normalized', 'VerticalAlignment', 'top', ...
  'FontName', 'Menlo', 'FontSize', 10, 'Interpreter', 'none');

exportgraphics(fig, pdfPath, 'ContentType', 'vector');
close(fig);
end

% -------------------------------------------------------------------------
function formatOffsetAxis(ax, offsets)
xlim(ax, [0.5 numel(offsets)+0.5]);
xticks(ax, 1:numel(offsets));
labels = arrayfun(@(x) sprintf('%g%c', x, char(176)), offsets, ...
  'UniformOutput', false);
xticklabels(ax, labels);
xlabel(ax, 'Probe offset');
end

% -------------------------------------------------------------------------
function jitter = deterministicJitter(n, width)
if n <= 1
  jitter = zeros(n,1);
else
  jitter = linspace(-width, width, n)';
end
end
