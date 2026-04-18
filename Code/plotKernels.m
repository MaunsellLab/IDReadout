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
%     kernels       : 5 x 2 x 2 x nVFrames kernels (diff/c/nC/RF/Opp, inc/dec steps, pref/probe)
%     kVars         : 5 x 2 x 2  (fixed over time) (diff/c/nC/RF/Opp, inc/dec steps, pref/probe)
%     compStats     : struct with 5 x 2 x 2 fields for kernel integrals and scale
%     nHits         : 5 x 2 hits, diff/c/nC/RF/Opp, dec/inc
%     nTrials       : 5 x 2 trials, diff/c/nC/RF/Opp, dec/inc
%
%   compStats.scale is the ordinary point estimate.
%   compStats.scaleCI, when present, provides bootstrap percentile intervals.
%
frameRateHz = header.frameRateHz.data;
stepMS = header.stepMS.data(1);
msPerVFrame = 1000.0 / frameRateHz;
[preStepMS, intStartMS, intDurMS] = integralWindowMS();
intStartVF = round((preStepMS + intStartMS) / msPerVFrame);
intStopVF = intStartVF + round(intDurMS / msPerVFrame);
trialDurMS = preStepMS + stepMS;
m = round(trialDurMS / msPerVFrame);
preStepVF = preStepMS / msPerVFrame;
plotTitles = {"Decrements", "Increments"};

% nSideTypes = size(kernels, 1);
nSideTypes = 5;   % don't plot the RF/Opp kernels
typeTitles = { ...
  '(Change RDK - No Change RDK)', 'Change RDK', 'No Change RDK', 'Left RDK', 'Right RDK', 'RF RDK', 'Opp RDK'};

figure(fig);
clf;
t = tiledlayout(nSideTypes, 2);
ax = nan(nSideTypes, 2);
overallYMax = 0;
overallYMin = 0;
prefLine = zeros(1, 2);
probeLine = zeros(1, 2);

for sideType = 1:nSideTypes
  for s = 1:2
    ax(sideType, s) = nexttile((sideType - 1) * 2 + 3 - s);
    xValues = 1:size(kernels, 4);
    hold on;
    [~, prefLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 1, :)), ...
      sqrt(kVars(sideType, s, 1)), [0, 0, 1]);
    [~, probeLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 2, :)), ...
      sqrt(kVars(sideType, s, 2)), [0.7, 0, 0.7]);
    plot([1, xValues(end)], [0, 0], 'k-');
    xticks([1, round(250 / msPerVFrame), round(500 / msPerVFrame), preStepVF, m]);
    xticklabels({"0", "250", "500", sprintf("%.0f", preStepMS), sprintf("%.0f", trialDurMS)});
    xlabel("Time (ms)");
    ylabel("Coherence (%)");
    title(sprintf("Coherence %s: %s Kernels", plotTitles{s}, typeTitles{sideType}));

    if isfield(compStats, 'scaleCI')
      textStr = sprintf(['Integrals: 0°: %.2f%%, ±%d°: %.2f%%\n' ...
        '\\fontsize{6}(integral: %.2f  scale = %.2f [SEM: %.2f; 68%% CI: %.2f, %.2f])'], ...
        compStats.kIntegrals(sideType, s, 1), probeDirDeg, compStats.kIntegrals(sideType, s, 2), ...
        compStats.kIntegrals(sideType, s, 2) / compStats.kIntegrals(sideType, s, 1), ...
        compStats.scale(sideType, s),compStats.scaleSEM(sideType, s), ...
        compStats.scaleCI.lo(sideType, s), ...
        compStats.scaleCI.hi(sideType, s));
    else
      textStr = sprintf(['Integrals: 0°: %.2f%%, ±%d°: %.2f%%\n' ...
        '\\fontsize{6}(integral: %.2f  scale = %.2f)'], ...
        compStats.kIntegrals(sideType, s, 1), probeDirDeg, compStats.kIntegrals(sideType, s, 2), ...
        compStats.kIntegrals(sideType, s, 2) / compStats.kIntegrals(sideType, s, 1), ...
        compStats.scale(sideType, s));
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
for sideType = 1:nSideTypes
  for a = 1:2
    ylim(ax(sideType, a), newYL);
    axes(ax(sideType, a)); %#ok<*LAXES>
    patch([preStepVF m m preStepVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    patch([intStartVF intStopVF intStopVF intStartVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
  end
end
axes(ax(1, 2));
legend([prefLine(1), probeLine(1)], {"0° ±SEM", sprintf("±%d°", probeDirDeg)}, 'location', 'southwest');
strElement = {"Decrements", "Increments"};
for s = 1:2
  subTitleStr{s} = sprintf('\\fontsize{8}%s: %d trials, %.0f%% correct, %d Left trials (%.1f%%), %.0f%% correct', ...
    strElement{3 - s}, hitStats.nTrials(3 - s), ...
    hitStats.nHits(3 - s) * 100.0 / hitStats.nTrials(3 - s), ...
    hitStats.nLeftTrials(3 - s), ...
    hitStats.nLeftTrials(3 - s) * 100 / hitStats.nTrials(3 - s), ...
    hitStats.nLeftHits(3 - s) * 100.0 / hitStats.nLeftTrials(3 - s)); %#ok<AGROW>
end
if ~isempty(probeDirDeg)
  probeStr = sprintf(' (Probe %d°)', probeDirDeg);
else
  probeStr = '';
end
title(t, sprintf('\\bf\\DeltaMean Kernels %s%s\\rm\n%s;         %s', ...
  underscoreToDash(titleStr), probeStr, subTitleStr{1}, subTitleStr{2}));
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

%-- underscoreToDash() --
function outStr = underscoreToDash(inStr)
% Replace underscores with dashes
if ischar(inStr)
  inStr = string(inStr);
end
outStr = replace(inStr, "_", "-");
end