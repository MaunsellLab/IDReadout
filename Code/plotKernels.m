function plotKernels(fig, titleStr, header, kernels, kVars, compStats, hitStats, probeDirDeg)
% Make a separate tile for inc/dec kernels, superimposing both preferred and
% probe kernel on each.

if nargin < 8
    probeDirDeg = [];
end
% Make a separate tile for inc/dec kernels, superimposing both preferred and
% probe kernel on each.
%
%   INPUTS:
%     fig           : number of the figure to use
%     header        : IDR header structure from dat file
%     kernels       : k x 2 x 2 x nVFrames kernels (numKTypes, dec/inc, pref/probe)
%     kVars         : k x 2 x 2  (fixed over time) (numKTypes, dec/inc, pref/probe)
%     compStats     : struct with 5 x 2 x 2 fields for kernel integrals and scale
%     nHits         : k x 2 hits, numKTypes, dec/inc
%     nTrials       : k x 2 trials, numKTypes, dec/inc
%     probeDirDeg   : probe direction for these kernels
%
%   Kernel traces are plotted as measured.
%   Text annotations use normalized comparison statistics when available:
%     compStats.normIntegrals, compStats.normR, compStats.normScale
%   with bootstrap intervals from compStats.normScaleCI when present.
%   Falls back to legacy raw fields for older compStats structs.
%

if nargin < 1 || isempty(fig)
    fig = 1;
end
if isnumeric(fig)   % Normalize fig to a figure handle
    figH = figure(fig);  % legacy behavior
elseif isgraphics(fig, 'figure')
    figH = fig;
else
    error('plotKernels:BadFigure: fig must be empty, a figure number, or a valid figure handle.');
end

frameRateHz = header.frameRateHz.data;
stepMS = header.stepMS.data(1);
msPerVFrame = 1000.0 / frameRateHz;
[preStepMS, intStartMS, intDurMS] = integralWindowMS();
intStartVF = round((preStepMS + intStartMS) / msPerVFrame);
intStopVF = intStartVF + round(intDurMS / msPerVFrame);
trialDurMS = preStepMS + stepMS;
trialDurVF = round(trialDurMS / msPerVFrame);
preStepVF = preStepMS / msPerVFrame;
plotTitles = {"Decrement", "Increment"};

nKernelTypes = size(kernels, 1);
typeTitles = {'(Change - No Change RDK)', 'Change RDK', 'No Change RDK', 'Left RDK', 'Right RDK', 'RF RDK', 'Opp RDK'};

set(figH, 'WindowStyle', 'docked');
clf(figH);
t = tiledlayout(figH, nKernelTypes, 2);
ax = nan(nKernelTypes, 2);
overallYMax = 0;
overallYMin = 0;
prefLine = zeros(1, 2);
probeLine = zeros(1, 2);
strElement = {"Dec", "Inc"};

for sideType = 1:nKernelTypes
  for s = 1:2
    colTitle = sprintf('%s: %d trials, %.0f%% correct (%d Left trials (%.1f%%), %.0f%% correct)', ...
      strElement{s}, hitStats.nTrials(3 - s), hitStats.nHits(3 - s) * 100.0 / hitStats.nTrials(3 - s), ...
      hitStats.nLeftTrials(3 - s), hitStats.nLeftTrials(3 - s) * 100 / hitStats.nTrials(3 - s), ...
      hitStats.nLeftHits(3 - s) * 100.0 / hitStats.nLeftTrials(3 - s));
    ax(sideType, s) = nexttile((sideType - 1) * 2 + 3 - s);
    xValues = 1:size(kernels, 4);
    hold on;
    [~, prefLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 1, :)), ...
      sqrt(kVars(sideType, s, 1)), [0, 0, 1]);
    [~, probeLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 2, :)), ...
      sqrt(kVars(sideType, s, 2)), [0.7, 0, 0.7]);
    plot([1, xValues(end)], [0, 0], 'k-');
    xticks([1, round(250 / msPerVFrame), round(500 / msPerVFrame), preStepVF, trialDurVF]);
    xticklabels({"0", "250", "500", sprintf("%.0f", preStepMS), sprintf("%.0f", trialDurMS)});
    xlabel("Time (ms)");
    ylabel("Coherence (%)");
    if sideType == 1
      title(sprintf("%s\n\n%s %s Kernels", colTitle, plotTitles{s}, typeTitles{sideType}));     
    else
      title(sprintf("%s %s Kernels", plotTitles{s}, typeTitles{sideType}));
    end
    % Text annotations report normalized comparison values when available.
    % The plotted kernel traces remain the raw measured kernels.
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
      % Backward compatibility for older compStats structs.
      displayIntegrals = compStats.kIntegrals;
      displayR         = compStats.R;
      displayScale     = compStats.scale;

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

    if ~isempty(displayScaleCI)
      textStr = sprintf(['\\fontsize{6}%s: 0°: %.2f%%, ±%d°: %.2f%%\n' ...
        '\\fontsize{6}(%s: %.2f; %s: %.2f [SEM: %.2f; 68%% CI: %.2f, %.2f])'], ...
        valueLabel, displayIntegrals(sideType, s, 1), probeDirDeg, displayIntegrals(sideType, s, 2), ...
        ratioLabel, displayR(sideType, s), ...
        scaleLabel, displayScale(sideType, s), displayScaleSEM(sideType, s), ...
        displayScaleCI.lo(sideType, s), ...
        displayScaleCI.hi(sideType, s));
    else
      textStr = sprintf(['\\fontsize{6}%s: 0°: %.2f%%, ±%d°: %.2f%%\n' ...
        '\\fontsize{6}(%s: %.2f; %s: %.2f)'], ...
        valueLabel, displayIntegrals(sideType, s, 1), probeDirDeg, displayIntegrals(sideType, s, 2), ...
        ratioLabel, displayR(sideType, s), ...
        scaleLabel, displayScale(sideType, s));
    end
    text(0.02, 0.98, textStr, 'units', 'normalized', ...
      'VerticalAlignment', 'top', 'fontSize', 9);
  end
  yl1 = ylim(ax(sideType, 1));
  yl2 = ylim(ax(sideType, 2));
  overallYMin = min(overallYMin, min(yl1(1), yl2(1)));
  overallYMax = max(overallYMax, max(yl1(2), yl2(2)));
end
set(ax, 'YLim', [overallYMin, overallYMax]);
newYL = [overallYMin, overallYMax];
for sideType = 1:nKernelTypes
  for a = 1:2
    ylim(ax(sideType, a), newYL);
    axes(ax(sideType, a)); %#ok<*LAXES>
    patch([preStepVF trialDurVF trialDurVF preStepVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    patch([intStartVF intStopVF intStopVF intStartVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
  end
end
axes(ax(1, 2));
legend([prefLine(1), probeLine(1)], {"0° ±SEM", "Probe"}, 'location', 'southwest');
title(t, underscoreToDash(titleStr));
end

%-- plotWithConstSEM() --
function [hPatch, hLine] = plotWithConstSEM(x, y, sem, faceColor)
  holdState = ishold;
  hold on;
  yUpper = y + sem;
  yLower = y - sem;
  xPatch = [x(:)', fliplr(x(:)')];
  yPatch = [yUpper(:)', fliplr(yLower(:)')];
  faceAlpha = 0.35;
  hPatch = patch(xPatch, yPatch, [0.7 0.7 0.7], 'faceColor', faceColor, 'FaceAlpha', faceAlpha, ...
    'EdgeColor', 'none');
  hLine = plot(x, y, '-' , 'color', faceColor, 'lineWidth', 1.5);
  uistack(hLine, 'top');
  
  if ~holdState
    hold off;
  end
end

%%
function outStr = underscoreToDash(inStr)
% Replace underscores with dashes
if ischar(inStr)
  inStr = string(inStr);
end
outStr = replace(inStr, "_", "-");
end

