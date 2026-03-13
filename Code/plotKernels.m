function plotKernels(fig, titleStr, header, kernels, kVars, kIntegrals, R, RVar, nHits, nTrials)
% Make a separate subplot for inc/dec kernels, superimposing both preferred and
% probe kernel on each. 
%
%   INPUTS:
%     fig           : number of the figure to use
%     header        : IDR header structure from dat file
%     kernels       : 5 x 2 x 2 x nVFrames kernels (diff/c/nC/RF/Opp, inc/dec steps, pref/probe)
%     kVars         : 5 x 2 x 2  (fixed over time) (diff/c/nC/RF/Opp, inc/dec steps, pref/probe) 
%     kIntegrals    : 5 x 2 x 2 kernel integrals during summing window
%     nHits         : 5 x 2 hits, diff/c/nC/RF/Opp, dec/inc
%     nTrials       : 5 x 2 trials, diff/c/nC/RF/Opp, dec/inc
%
  frameRateHz = header.frameRateHz.data;
  stepMS = header.stepMS.data(1);
  probeDirDeg = header.probeDirDeg.data(1);
  msPerVFrame = 1000.0 / frameRateHz;
  [preStepMS, intStartMS, intDurMS] = integralWindowMS();
  intStartVF = round((preStepMS + intStartMS) / msPerVFrame);
  intStopVF = intStartVF + round(intDurMS / msPerVFrame);
  trialDurMS = preStepMS + stepMS;
  m = round(trialDurMS / msPerVFrame); 
  preStepVF = preStepMS / msPerVFrame;
  plotTitles = {"Decrements", "Increments"};
  typeTitles = {"Difference", "Change Side", "No Change Side", "RF Side", "Opp Side"};

  figure(1);
  clf;
  ax = nan(5, 2);
  overallYMax = 0;
  overallYMin = 0;
  prefLine = zeros(1, 2);
  probeLine = zeros(1, 2);
  for sideType = 1:5
    for s = 1:2                           % for each coherence step direction (inc/dec)    
      ax(sideType, s) = subplot(5, 2, (sideType - 1) * 2 + 3 - s);
      xValues = 1:size(kernels, 4);
      hold on;
      [~, prefLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 1, :)), ...
                sqrt(kVars(sideType, s, 1)), [0, 0, 1]);
      [~, probeLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(sideType, s, 2, :)), ...
                sqrt(kVars(sideType, s, 2)), [0.7, 0, 0.7]);
      plot([1, xValues(end)], [0, 0], 'k-');
      xticks([1, round(250 / msPerVFrame), round(500 / msPerVFrame), preStepVF, m]);
      xticklabels({"0", "250", "500", sprintf("%.0f", preStepMS), sprintf("%.0f", trialDurMS)});
      xlabel("Time During Trial (ms)");
      ylabel("Percent Coherence");
      title(sprintf("%s: \\DeltaMean %s %s Kernels", titleStr, typeTitles{sideType}, plotTitles{s}));
      text(0.02, 0.98, sprintf("0°: %.2f%%, ±%d°: %.2f%% (R = %.2f ± %.2f SEM)\nn = %d (%.0f%% correct)", ...
          kIntegrals(sideType, s, 1), probeDirDeg, kIntegrals(sideType, s, 2), R(sideType, s), ...
                sqrt(RVar(sideType, s)), nTrials(sideType, s), ...
          nHits(sideType, s) * 100.0 / nTrials(sideType, s)), ...
          'units', 'normalized', 'VerticalAlignment', 'top', 'fontSize', 9);
    end
    yl1 = ylim(ax(sideType, 1));
    yl2 = ylim(ax(sideType, 2));
    overallYMin = min(overallYMin, min(yl1(1), yl2(1)));
    overallYMax = max(overallYMax, max(yl1(2), yl2(2)));
  end
  set(ax, 'YLim', [overallYMin, overallYMax]);
  newYL = [overallYMin, overallYMax];
  for sideType = 1:5
    for a = 1:2
      ylim(ax(sideType, a), newYL);
      axes(ax(sideType, a));
      patch([preStepVF m m preStepVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
      patch([intStartVF intStopVF intStopVF intStartVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
        [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
    end
  end
  axes(ax(1, 2));
  legend([prefLine(1), probeLine(1)], {"0° ±SEM", sprintf("±%d°", probeDirDeg)}, 'location', 'southwest');


end

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

