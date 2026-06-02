function plotSideTypeKernelAverage(varargin)
% plotSideTypeKernelAverage
%
% Plot average kernels for one side type across probe directions.
%
% This is a downstream plotting function.  It does not recompute kernels or
% averages.  It loads AverageKernelPlotData.mat files saved by kernelAverage().
%
% Example:
%   plotSideTypeKernelAverage()
%   plotSideTypeKernelAverage('SideType', 'change')
%   plotSideTypeKernelAverage('SideType', 'noChange')
%   plotSideTypeKernelAverage('ProbeDirs', [45 90 135 180])

baseFolder = folderPath();

P = inputParser;
addParameter(P, 'SideType', 'change', @(x) ischar(x) || isstring(x));
addParameter(P, 'ProbeDirs', [], @(x) isempty(x) || isnumeric(x));
parse(P, varargin{:});
R = P.Results;

sideType = char(R.SideType);
sideTypeIndex(sideType);  % validate sideType name; index not otherwise needed here

if isempty(R.ProbeDirs)
  probeDirs = findProbeDirsWithAverageData(sideType);
else
  probeDirs = R.ProbeDirs(:).';
end

assert(~isempty(probeDirs), 'plotSideTypeKernelAverage:NoProbeDirs', ...
  'No probe directions found for sideType "%s".', sideType);

% ---- Load average-kernel plot data ----
nProbes = numel(probeDirs);
plotData = cell(1, nProbes);
missingIndices = [];

for p = 1:nProbes
  probeTag = sprintf('Probe%d', probeDirs(p));
  dataFile = fullfile(baseFolder, 'Data', probeTag, 'AverageKernels', sideType, ...
          'AverageKernelPlotData.mat');
  if ~isfile(dataFile)
    missingIndices = [missingIndices, p]; %#ok<AGROW>
    fprintf('      no kernels for probe %3d°\n', probeDirs(p));
    continue;
  end

  S = load(dataFile);
  D = S.averageKernelPlotData;
  assert(isfield(D, 'summarySideType'), 'plotSideTypeKernelAverage:MissingSideType', ...
    'averageKernelPlotData.summarySideType missing in %s', dataFile);
  assert(strcmp(D.summarySideType, sideType), 'plotSideTypeKernelAverage:WrongSideType', ...
    'Requested sideType "%s", but %s contains sideType "%s".', ...
    sideType, dataFile, D.summarySideType);
  assert(isfield(D, 'kernels') && ndims(D.kernels) == 3, 'plotSideTypeKernelAverage:BadKernels', ...
    'Expected kernels to be step x stream x time in %s', dataFile);
  assert(isfield(D, 'stepTypeNames') && numel(D.stepTypeNames) == size(D.kernels, 1), ...
    'plotSideTypeKernelAverage:BadStepTypeNames', 'stepTypeNames missing or inconsistent in %s', dataFile);
  assert(isfield(D, 'streamTypeNames') && numel(D.streamTypeNames) == size(D.kernels, 2), ...
    'plotSideTypeKernelAverage:BadStreamTypeNames', 'streamTypeNames missing or inconsistent in %s', dataFile);
  assert(isfield(D, 'tMS') && numel(D.tMS) == size(D.kernels, 3), 'plotSideTypeKernelAverage:BadTimeVector', ...
    'tMS missing or inconsistent in %s', dataFile);
  plotData{p} = D;
end
if ~isempty(missingIndices)
  probeDirs(missingIndices) = [];
  plotData(missingIndices) = [];
  nProbes = numel(probeDirs);
end

% ---- Confirm compatible labels and timing ----
refStepTypeNames = plotData{1}.stepTypeNames;
refStreamTypeNames = plotData{1}.streamTypeNames;
refTMS = plotData{1}.tMS(:).';

for p = 2:nProbes
  assert(isequal(plotData{p}.stepTypeNames, refStepTypeNames), ...
    'plotSideTypeKernelAverage:StepTypeMismatch', 'stepTypeNames differ across probe directions.');
  assert(isequal(plotData{p}.streamTypeNames, refStreamTypeNames), ...
    'plotSideTypeKernelAverage:StreamTypeMismatch', 'streamTypeNames differ across probe directions.');
  assert(numel(plotData{p}.tMS) == numel(refTMS) &&  all(abs(plotData{p}.tMS(:).' - refTMS) < 1e-9), ...
    'plotSideTypeKernelAverage:TimeMismatch', 'tMS differs across probe directions.');
end
nSteps = numel(refStepTypeNames);
nStreams = numel(refStreamTypeNames);

% ---- Common y-limits, including SEM patches ----
overallYMin = zeros(1, nProbes);
overallYMax = zeros(1, nProbes);

for p = 1:nProbes
  D = plotData{p};
  for stepType = 1:nSteps
    for streamType = 1:nStreams
      y = squeeze(D.kernels(stepType, streamType, :));

      if isfield(D, 'kVars') && ~isempty(D.kVars)
        sem = sqrt(D.kVars(stepType, streamType));
      else
        sem = 0;
      end

      yUpper = y + sem;
      yLower = y - sem;

      finiteVals = [yUpper(:); yLower(:)];
      finiteVals = finiteVals(isfinite(finiteVals));

      if ~isempty(finiteVals)
        overallYMin(p) = min(overallYMin(p), min(finiteVals));
        overallYMax(p) = max(overallYMax(p), max(finiteVals));
      end
    end
  end
end

% ---- Timing for x-axis and gray patches ----
  D0 = plotData{1};  
  msPerVFrame = D0.msPerVFrame;
  preStepMS   = D0.firstPreStepMS;
  stepMS      = D0.firstStepMS;
  [~, intStartMS, intDurMS] = integralWindowMS();
  
  intStartVF = round((preStepMS + intStartMS) / msPerVFrame);
  intStopVF  = intStartVF + round(intDurMS / msPerVFrame);
  
  trialDurMS = preStepMS + stepMS;
  trialDurVF = round(trialDurMS / msPerVFrame);
  preStepVF  = preStepMS / msPerVFrame;
  
  nFrames = size(plotData{1}.kernels, 3);
  xValues = 1:nFrames;
  
  xTickVals = [1, round(250 / msPerVFrame), round(500 / msPerVFrame), ...
    preStepVF, trialDurVF];
  
  xTickLabels = {"0", "250", "500", ...
    sprintf("%.0f", preStepMS), sprintf("%.0f", trialDurMS)};

% ---- Plot ----
fig = figure(200);
t = tiledlayout(fig, nProbes, nSteps, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(nProbes, nSteps);
prefLine = gobjects(1, nSteps);
probeLine = gobjects(1, nSteps);
plotTitles = {'Decrement', 'Increment'};

for p = 1:nProbes
  newYL = [overallYMin(p) overallYMax(p)];
  D = plotData{p};
  for s = 1:nSteps
    % Match plotKernels convention: dec/inc are indexed 1/2 but displayed
    % with increment on the left and decrement on the right.  We flip the
    % column locations here, but everything else should stay straight
    % indexing
    ax(p, s) = nexttile((p - 1) * nSteps + 3 - s);
    hold on;
    [~, prefLine(s)] = plotWithConstSEM(xValues, squeeze(D.kernels(s, 1, :)), sqrt(D.kVars(s, 1)), [0, 0, 1]);
    [~, probeLine(s)] = plotWithConstSEM(xValues, squeeze(D.kernels(s, 2, :)), sqrt(D.kVars(s, 2)), [0.7, 0, 0.7]);
    plot([1, xValues(end)], [0, 0], 'k-');
    plot([1, xValues(end)], [0, 0], 'k-');
    yline(1, ':', 'Color', [0.5 0.5 0.5]);
    yline(-1, ':', 'Color', [0.5 0.5 0.5]);
    % Gray patches.
    patch([preStepVF trialDurVF trialDurVF preStepVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    patch([intStartVF intStopVF intStopVF intStartVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    textStr = makeKernelStatsText(D, s);
    text(0.02, 0.98, textStr, 'units', 'normalized',  'VerticalAlignment', 'top',  'fontSize', 6);

    ylim(newYL);
    xlim([1 xValues(end)]);
    xticks(xTickVals);
    xticklabels(xTickLabels);
    xlabel("Time (ms)");
    ylabel("Coherence (%)");
    title(sprintf('%s Probe %s. (%d Sessions, %d Trials)', D.titlePrefix, plotTitles{s}(1:3), D.nSessions, ...
      D.avgHitStats.nTrials(s)));
    box on;
  end
end

upperSideType = [upper(sideType(1)) sideType(2:end)];
legend([prefLine(1), probeLine(1)], {"0° ±SEM", "Probe"}, 'location', 'southwest');
title(t, sprintf('\\bf%s Side Average Kernels by Probe Direction', upperSideType));

% ---- Export ----
outFolder = fullfile(baseFolder, 'Plots', 'AverageKernels', upperSideType);
validFolder(outFolder);

outFile = fullfile(outFolder, sprintf('SideTypeKernelAverage_%s.pdf', sideType));
exportgraphics(fig, outFile, 'ContentType', 'vector');
end

%% Find probe directions with saved average-kernel data for requested side type.
function probeDirs = findProbeDirsWithAverageData(sideType)

dataFolder = fullfile(folderPath(), 'Data');
probeFolders = dir(fullfile(dataFolder, 'Probe*'));

probeDirs = [];

for i = 1:numel(probeFolders)
  if ~probeFolders(i).isdir
    continue;
  end

  probeName = probeFolders(i).name;
  token = regexp(probeName, '^probe(-?\d+)$', 'tokens', 'once');

  if isempty(token)
    continue;
  end

  probeDir = str2double(token{1});
  dataFile = fullfile(dataFolder, probeName, 'AverageKernels', sideType, 'AverageKernelPlotData.mat');

  if isfile(dataFile)
    probeDirs(end+1) = probeDir; %#ok<AGROW>
  end
end

probeDirs = sort(probeDirs);

end

%% Make statistics for one kernel plots
function textStr = makeKernelStatsText(D, stepType)

compStats = D.avgCompStats;
sideTypeNum = D.sideTypeNum;
probeDirDeg = D.probeDirDeg;

% Prefer normalized comparison statistics. The plotted traces remain the
% raw measured kernels, but the displayed integral, ratio, and scale values
% are the amplitude-normalized values used for interpretation.
if isfield(compStats, 'normIntegrals')
  displayIntegrals = compStats.normIntegrals;
  displayR         = compStats.normR;
  displayScale     = compStats.normScale;

  if isfield(compStats, 'normScaleSEM')
    displayScaleSEM = compStats.normScaleSEM;
  else
    displayScaleSEM = nan(size(displayScale));
  end

  if isfield(compStats, 'normScaleCI')
    displayScaleCI = compStats.normScaleCI;
  else
    displayScaleCI = [];
  end

  valueLabel = 'norm integrals';
  ratioLabel = 'norm ratio';
  scaleLabel = 'norm scale';
else
  % Backward compatibility for older AverageKernelPlotData files.
  displayIntegrals = compStats.kIntegrals;

  if isfield(compStats, 'R')
    displayR = compStats.R;
  else
    displayR = nan(size(compStats.scale));
    for iSide = 1:size(displayIntegrals, 1)
      for iStep = 1:size(displayIntegrals, 2)
        prefIntegral = displayIntegrals(iSide, iStep, 1);
        probeIntegral = displayIntegrals(iSide, iStep, 2);
        if abs(prefIntegral) > eps
          displayR(iSide, iStep) = probeIntegral / prefIntegral;
        end
      end
    end
  end

  displayScale = compStats.scale;

  if isfield(compStats, 'scaleSEM')
    displayScaleSEM = compStats.scaleSEM;
  else
    displayScaleSEM = nan(size(displayScale));
  end

  if isfield(compStats, 'scaleCI')
    displayScaleCI = compStats.scaleCI;
  else
    displayScaleCI = [];
  end

  valueLabel = 'integrals';
  ratioLabel = 'ratio';
  scaleLabel = 'scale';
end

prefIntegral  = displayIntegrals(sideTypeNum, stepType, 1);
probeIntegral = displayIntegrals(sideTypeNum, stepType, 2);
integralRatio = displayR(sideTypeNum, stepType);
scaleVal      = displayScale(sideTypeNum, stepType);
scaleSEM      = displayScaleSEM(sideTypeNum, stepType);

if ~isempty(displayScaleCI)
  textStr = sprintf(['%s: 0°: %.2f%%, %s: %.2f%%\n' '(%s: %.2f; %s: %.2f [SEM: %.2f; 68%% CI: %.2f, %.2f])'], ...
    valueLabel, prefIntegral, probeDirString(probeDirDeg), probeIntegral, ratioLabel, integralRatio, scaleLabel, ...
    scaleVal, scaleSEM, displayScaleCI.lo(sideTypeNum, stepType), displayScaleCI.hi(sideTypeNum, stepType));
else
  textStr = sprintf(['%s: 0°: %.2f%%, %s: %.2f%%\n' '%s: %.2f; %s: %.2f)'], ...
    valueLabel, prefIntegral, probeDirString(probeDirDeg), probeIntegral, ratioLabel, integralRatio, ...
    scaleLabel, scaleVal);
end
end

%% probeDirString() -- return a string with appropriate ± and degree symbol
function probeStr = probeDirString(probeDir)

  if (probeDir == 180 || probeDir == 0) 
    probeStr = sprintf('%d°', probeDir);
  else
    probeStr = sprintf('±%d°', probeDir);
  end
end

%% plotWithConstSEM
function [hPatch, hLine] = plotWithConstSEM(x, y, sem, faceColor)

holdState = ishold;
  hold on;
  yUpper = y + sem;
  yLower = y - sem;
  xPatch = [x(:)', fliplr(x(:)')];
  yPatch = [yUpper(:)', fliplr(yLower(:)')];
  faceAlpha = 0.35;
  hPatch = patch(xPatch, yPatch, [0.7 0.7 0.7], 'faceColor', faceColor, 'FaceAlpha', faceAlpha, 'EdgeColor', 'none');
  hLine = plot(x, y, '-', 'color', faceColor, 'lineWidth', 1.5);
  uistack(hLine, 'top');
  if ~holdState
    hold off;
  end
end