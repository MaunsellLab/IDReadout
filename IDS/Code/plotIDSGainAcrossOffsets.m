function [gainTable, fig] = plotIDSGainAcrossOffsets(summaries)
% plotIDSGainAcrossOffsets  Plot IDS probe gains across tested offsets.
%
%   gainTable = plotIDSGainAcrossOffsets(summaries)
%   [gainTable, fig] = plotIDSGainAcrossOffsets(summaries)
%
% summaries is the output from plotIDSSummary. The function saves an
% across-offset table (MAT and CSV) and, when at least two offsets are
% present, creates a docked/display figure and a portrait PDF export. The
% PDF name includes the animal selection stored by plotIDSSummary, for
% example IDSGainAcrossOffsets_Meetz.pdf or IDSGainAcrossOffsets_All.pdf.

sigmaMTDeg = 37.5;
root = domainFolder(mfilename('fullpath'));
summaryFolder = validFolder(fullfile(root, 'Data', 'AcrossOffsetSummaries'));
plotFolder = validFolder(fullfile(root, 'Plots'));

gainTable = makeGainTable(summaries);
fig = gobjects(0);

if isempty(gainTable)
  warning('plotIDSGainAcrossOffsets:NoSummaries', ...
    'No IDS offset summaries were supplied; no files were created.');
  return
end

animal = summaryAnimal(summaries);

metadata = struct();
metadata.createdBy = mfilename;
metadata.createdAt = datetime('now');
metadata.animal = animal;
metadata.mtReferenceSigmaDeg = sigmaMTDeg;
metadata.absoluteMTReferenceScale = mean(gainTable.driftGain, 'omitnan');
metadata.absoluteMTReferenceScaleMethod = ...
  'Unweighted mean drift gain across probe offsets';
metadata.normalizedGain = 'Pooled probe gain divided by drift gain';
metadata.ratioCIMethod = char(gainTable.ratioCIMethod(1));

save(fullfile(summaryFolder, 'IDS_AcrossOffsetGainSummary.mat'), ...
  'gainTable', 'metadata');
writetable(gainTable, ...
  fullfile(summaryFolder, 'IDS_AcrossOffsetGainSummary.csv'));

if height(gainTable) < 2
  warning('plotIDSGainAcrossOffsets:OnlyOneOffset', ...
    ['Saved the across-offset gain table, but at least two offsets are ' ...
     'required to create the across-offset figure.']);
  return
end

% The ordinary figure respects the user's default docked-window setting.
fig = figure('Color', 'w');
plotGainFigure(fig, gainTable, sigmaMTDeg);

% Render independently into a fixed portrait canvas for PDF export.
exportFig = figure('Color', 'w', 'WindowStyle', 'normal', 'Visible', 'off', 'Units', 'inches', ...
  'Position', [1 1 8 10.5]);
cleanup = onCleanup(@() deleteFigure(exportFig)); 
plotGainFigure(exportFig, gainTable, sigmaMTDeg);
drawnow;

pdfName = sprintf('IDSGainAcrossOffsets_%s.pdf', animal);
exportgraphics(exportFig, fullfile(plotFolder, pdfName), ...
  'ContentType', 'vector', 'Units', 'inches', 'Width', 8, 'Height', 10.5, 'Padding', 0.1, 'PreserveAspectRatio', 'on');
end

function animal = summaryAnimal(summaries)
if isstruct(summaries)
  entries = num2cell(summaries(:));
else
  entries = summaries(:);
end

if ~all(cellfun(@(S) isfield(S, 'animal'), entries))
  error('plotIDSGainAcrossOffsets:MissingAnimal', ...
    ['Each summary must contain the animal selection recorded by ' ...
     'plotIDSSummary. Recreate legacy summaries with plotIDSSummary.']);
end

animals = string(cellfun(@(S) S.animal, entries, 'UniformOutput', false));
animals = unique(animals, 'stable');
if numel(animals) ~= 1
  error('plotIDSGainAcrossOffsets:MixedAnimals', ...
    'All offset summaries must have the same animal selection.');
end

switch lower(animals)
  case "all"
    animal = 'All';
  case "meetz"
    animal = 'Meetz';
  case "neesha"
    animal = 'Neesha';
  otherwise
    error('plotIDSGainAcrossOffsets:UnknownAnimal', ...
      'Unknown animal selection %s.', animals);
end
end

function T = makeGainTable(summaries)
if isempty(summaries)
  T = table();
  return
end

if isstruct(summaries)
  entries = num2cell(summaries(:));
elseif iscell(summaries) && all(cellfun(@isstruct, summaries(:)))
  entries = summaries(:);
else
  error('plotIDSGainAcrossOffsets:BadSummaries', ...
    'summaries must be a structure or cell array of structures.');
end

n = numel(entries);
probeOffsetDeg = nan(n,1);
nSessions = nan(n,1);
nTrials = nan(n,1);
nTrialsMeetz = nan(n,1);
nTrialsNeesha = nan(n,1);
nStepNoiseTrials = nan(n,1);
nStepNoiseTrialsMeetz = nan(n,1);
nStepNoiseTrialsNeesha = nan(n,1);
noisePredictorType = strings(n,1);
driftGain = nan(n,1);
driftGainCILow = nan(n,1);
driftGainCIHigh = nan(n,1);
pooledProbeGain = nan(n,1);
pooledProbeGainCILow = nan(n,1);
pooledProbeGainCIHigh = nan(n,1);
probeOverDrift = nan(n,1);
ratioCI95Low = nan(n,1);
ratioCI95High = nan(n,1);
ratioCIMethod = strings(n,1);

for k = 1:n
  S = entries{k};
  validateSummary(S, k);
  F = S.noiseGain.driftProbe;
  N = S.noiseGain.normalized;

  probeOffsetDeg(k) = S.probeOffsetDeg;
  nSessions(k) = S.nSessions;
  nTrials(k) = S.nTrials;
  nStepNoiseTrials(k) = S.nStepNoiseTrials;
  [nTrialsMeetz(k), nTrialsNeesha(k), nStepNoiseTrialsMeetz(k), nStepNoiseTrialsNeesha(k)] = animalTrialCounts(S);
  noisePredictorType(k) = string(S.noisePredictorType);
  driftGain(k) = F.gain(1);
  driftGainCILow(k) = F.CI95(1,1);
  driftGainCIHigh(k) = F.CI95(1,2);
  pooledProbeGain(k) = F.gain(2);
  pooledProbeGainCILow(k) = F.CI95(2,1);
  pooledProbeGainCIHigh(k) = F.CI95(2,2);
  probeOverDrift(k) = N.pooledProbeOverDrift;
  ratioCI95Low(k) = N.pooledCI95(1);
  ratioCI95High(k) = N.pooledCI95(2);
  ratioCIMethod(k) = string(N.method);
end

if numel(unique(noisePredictorType)) ~= 1
  error('plotIDSGainAcrossOffsets:MixedPredictors', ...
    'All offset summaries must use the same noise predictor type.');
end
if numel(unique(probeOffsetDeg)) ~= n
  error('plotIDSGainAcrossOffsets:DuplicateOffsets', ...
    'More than one summary was supplied for the same probe offset.');
end
if numel(unique(ratioCIMethod)) ~= 1
  error('plotIDSGainAcrossOffsets:MixedCIMethods', ...
    'All offset summaries must use the same normalized-gain CI method.');
end

T = table(probeOffsetDeg, nSessions, nTrials, nTrialsMeetz, nTrialsNeesha, ...
  nStepNoiseTrials, nStepNoiseTrialsMeetz, nStepNoiseTrialsNeesha, ...
  noisePredictorType, driftGain, driftGainCILow, driftGainCIHigh, ...
  pooledProbeGain, pooledProbeGainCILow, pooledProbeGainCIHigh, ...
  probeOverDrift, ratioCI95Low, ratioCI95High, ratioCIMethod);
T = sortrows(T, 'probeOffsetDeg');
end

function [nMeetz, nNeesha, nStepNoiseMeetz, nStepNoiseNeesha] = animalTrialCounts(S)
nMeetz = 0;
nNeesha = 0;
if ~isfield(S, 'sessionRecords') || ~istable(S.sessionRecords) || ...
    ~ismember('nTrials', S.sessionRecords.Properties.VariableNames) || ...
    ~ismember('nStepNoiseTrials', S.sessionRecords.Properties.VariableNames)
  return
end

R = S.sessionRecords;
if ismember('animal', R.Properties.VariableNames)
  animal = string(R.animal);
elseif ismember('sourceSessionID', R.Properties.VariableNames)
  sourceID = string(R.sourceSessionID);
  animal = strings(height(R), 1);
  animal(contains(sourceID, 'Meetz', 'IgnoreCase', true)) = "Meetz";
  animal(contains(sourceID, 'Neesha', 'IgnoreCase', true)) = "Neesha";
else
  return
end
nMeetz = sum(R.nTrials(strcmpi(animal, 'Meetz')), 'omitnan');
nNeesha = sum(R.nTrials(strcmpi(animal, 'Neesha')), 'omitnan');
nStepNoiseMeetz = sum(R.nStepNoiseTrials(strcmpi(animal, 'Meetz')), 'omitnan');
nStepNoiseNeesha = sum(R.nStepNoiseTrials(strcmpi(animal, 'Neesha')), 'omitnan');
end

function validateSummary(S, index)
required = {'probeOffsetDeg','nSessions','nTrials','nStepNoiseTrials', 'noisePredictorType','noiseGain'};
missing = required(~isfield(S, required));
if ~isempty(missing)
  error('plotIDSGainAcrossOffsets:IncompleteSummary', ...
    'Summary %d is missing field %s.', index, missing{1});
end
if ~isfield(S.noiseGain, 'driftProbe') || ...
    ~isfield(S.noiseGain, 'normalized')
  error('plotIDSGainAcrossOffsets:IncompleteGainFit', ...
    'Summary %d lacks the required driftProbe or normalized gain fit.', index);
end
F = S.noiseGain.driftProbe;
N = S.noiseGain.normalized;
if numel(F.gain) ~= 2 || ~isequal(size(F.CI95), [2 2]) || ...
    ~isfield(N, 'pooledProbeOverDrift') || ...
    ~isfield(N, 'pooledCI95') || numel(N.pooledCI95) ~= 2 || ...
    ~isfield(N, 'method')
  error('plotIDSGainAcrossOffsets:UnexpectedGainFit', ...
    'Summary %d has an unexpected gain-fit structure.', index);
end
end

%%-----------------------------------------
function plotGainFigure(fig, T, sigmaMTDeg)
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('IDS Probe Gain Across Offsets (%s predictor)', char(T.noisePredictorType(1))), 'FontWeight', 'bold');

offset = T.probeOffsetDeg;
xLimits = offsetLimits(offset);
xTicks = unique([0; offset; 180]);

ax = nexttile(tl, 1);
hold(ax, 'on');
xReference = linspace(0, 180, 721);
yReference = exp(-0.5 * (xReference / sigmaMTDeg).^2);
plot(ax, xReference, yReference, '--', 'Color', [.35 .35 .35], 'LineWidth', 1.5, ...
  'DisplayName', sprintf('MT tuning (\\sigma = %.1f^\\circ)', sigmaMTDeg));
plotAsymmetricErrorBars(ax, offset, T.probeOverDrift, T.ratioCI95Low, T.ratioCI95High, 'o-', [0.55 0.20 0.75], ...
  'Pooled probe/drift gain');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlim(ax, xLimits);
ylim(ax, [-0.5, 2.5]);
set(ax, 'XTick', xTicks);
xlabel(ax, 'Probe Offset (deg)');
ylabel(ax, 'Normalized Probe Gain  g_P/g_D');
title(ax, 'Drift-Normalized Probe Gain');
grid(ax, 'on');
box(ax, 'off');
legend(ax, 'Location', 'northwest', 'FontSize', 8);
addSummaryBox(ax, T);

ax = nexttile(tl, 2);
hold(ax, 'on');
absoluteMTScale = mean(T.driftGain, 'omitnan');
plot(ax, xReference, absoluteMTScale * yReference, '--', ...
  'Color', [.35 .35 .35], 'LineWidth', 1.5, ...
  'DisplayName', sprintf('MT tuning, mean g_D (%.3f)', ...
  absoluteMTScale));
plotAsymmetricErrorBars(ax, offset, T.driftGain, ...
  T.driftGainCILow, T.driftGainCIHigh, 'o-', [0.10 0.45 0.80], ...
  'Drift gain');
plotAsymmetricErrorBars(ax, offset, T.pooledProbeGain, ...
  T.pooledProbeGainCILow, T.pooledProbeGainCIHigh, 's-', [0.75 0.35 0.70], ...
  'Pooled probe gain');
yline(ax, 0, 'k:', 'HandleVisibility', 'off');
xlim(ax, xLimits);
ylim(ax, [-0.5, 2.5]);
set(ax, 'XTick', xTicks);
xlabel(ax, 'Probe Offset (deg)');
ylabel(ax, 'Absolute Noise Gain');
title(ax, 'Absolute Drift and Probe Gains');
grid(ax, 'on');
box(ax, 'off');
legend(ax, 'Location', 'northwest', 'FontSize', 8);
addAbsoluteSummaryBox(ax, T);
end

%%-----------------------------------------
function plotAsymmetricErrorBars(ax, x, y, low, high, lineSpec, color, label)
lowerLength = y - low;
upperLength = high - y;
errorbar(ax, x, y, lowerLength, upperLength, lineSpec, ...
  'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 5, ...
  'LineWidth', 1.3, 'CapSize', 8, 'DisplayName', label);
end

%%-----------------------------------------
function addSummaryBox(ax, T)
lines = strings(height(T) + 1, 1);
lines(1) = "Offset   g_P/g_D     95% CI    Meetz Neesha";
for k = 1:height(T)
  lines(k+1) = sprintf('%5.0f%c %6.2f [%5.2f, %5.2f] %5d %5d', ...
    T.probeOffsetDeg(k), char(176), T.probeOverDrift(k), ...
    T.ratioCI95Low(k), T.ratioCI95High(k), ...
    T.nStepNoiseTrialsMeetz(k), T.nStepNoiseTrialsNeesha(k));
end
text(ax, .98, .98, cellstr(lines), 'Units', 'normalized', ...
  'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
  'FontName', 'Menlo', 'FontSize', 7.5, 'BackgroundColor', 'w', ...
  'EdgeColor', [.65 .65 .65], 'Margin', 5);
end

%%-----------------------------------------
function addAbsoluteSummaryBox(ax, T)
lines = strings(height(T) + 1, 1);
lines(1) = "Offset   g_D    (95% CI)      g_P     (95% CI)";
for k = 1:height(T)
  lines(k+1) = sprintf('%5.0f%c %5.2f [%4.2f, %4.2f], %5.2f [%5.2f, %5.2f]', ...
    T.probeOffsetDeg(k), char(176), T.driftGain(k), T.driftGainCILow(k), T.driftGainCIHigh(k), ...
    T.pooledProbeGain(k), T.pooledProbeGainCILow(k), T.pooledProbeGainCIHigh(k));
end
text(ax, .98, .98, cellstr(lines), 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
  'VerticalAlignment', 'top', 'FontName', 'Menlo', 'FontSize', 7.5, 'BackgroundColor', 'w', ...
  'EdgeColor', [.65 .65 .65], 'Margin', 5);
end

%%-----------------------------------------
function limits = offsetLimits(offset)
maxOffset = max(180, max(offset));
limits = [0 maxOffset];
if all(offset == 0)
  limits = [0 1];
end
end

%%-----------------------------------------
function deleteFigure(fig)
if isgraphics(fig)
  delete(fig);
end
end
