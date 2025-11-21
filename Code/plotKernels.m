function plotKernels(fig, titleStr, header, kernels, kVars, kIntegrals, R, RVar, nHits, nTrials)
% Make a separate subplot for inc/dec kernels, superimposing both preferred and
% probe kernel on each. 
%
%   INPUTS:
%     fig           : number of the figure to use
%     header        : IDR header structure from dat file
%     kernels       : 2 x 2 x nVFrames kernel values (inc/dec steps,
%                     pref/probe, kernel values)
%     kVars         : 2 x 2 kernel variances (fixed over time) (inc/dec steps, pref/probe) 
%     kIntegrals    : 2 x 2 integrals of kernels during the summing window
%     intStartMS    : start of summing window
%     intStopMS     : end of summing widnow
%     nHits         : 1 x 2 vector of hits, dec/inc
%     nTrials       : 1 x 2 vector of trials, dec/inc
%
%   OUTPUTS:
%     kernel   : m x 1 kernel, defined as
%                  mean(noise | wrong) - mean(noise | correct)
%     kernelSE : m x 1 standard error of the kernel at each time bin (based on known distribution variance)
%
  frameRateHz = header.frameRateHz.data;
  % preStepMS = header.preStepMS.data(1);       % sometimes header values are replicated vectors
  stepMS = header.stepMS.data(1);
  probeDirDeg = header.probeDirDeg.data(1);
  msPerVFrame = 1000.0 / frameRateHz;
  [preStepMS, intStartMS, intDurMS] = integralWindowMS();
  intStartVF = round((preStepMS + intStartMS) / msPerVFrame);
  intStopVF = intStartVF + round(intDurMS / msPerVFrame);
  trialDurMS = preStepMS + stepMS;
  m = round(trialDurMS / msPerVFrame); 
  preStepVF = preStepMS / msPerVFrame;
  plotTitles = {"Decrement", "Increment"};

  figure(fig);
  clf
  ax = nan(1, 2);
  prefLine = zeros(1, 2);
  probeLine = zeros(1, 2);
  for s = 1:2                           % for each coherence step direction (inc/dec)    
    ax(s) = subplot(2, 1, 3 - s);
    xValues = 1:size(kernels, 3);
    hold on;
    [~, prefLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(s, 1, :)), sqrt(kVars(s, 1)), [0, 0, 1]);
    [~, probeLine(s)] = plotWithConstSEM(xValues, squeeze(kernels(s, 2, :)), sqrt(kVars(s, 2)), [0.7, 0, 0.7]);
    plot([1, xValues(end)], [0, 0], 'k-');
    xticks([1, round(250 / msPerVFrame), round(500 / msPerVFrame), preStepVF, m]);
    xticklabels({"0", "250", "500", sprintf("%.0f", preStepMS), sprintf("%.0f", trialDurMS)});
    xlabel("Time During Trial (ms)");
    ylabel("Percent Coherence");
    title(sprintf("%s: Coherence %s, \\DeltaMean Kernels", titleStr, plotTitles{s}));
    text(0.02, 0.98, sprintf("0°: %.2f%%, ±%d°: %.2f%% (R = %.2f ± %.2f SEM)\nn = %d (%.0f%% correct)", ...
        kIntegrals(s, 1), probeDirDeg, kIntegrals(s, 2), R(s), sqrt(RVar(s)), nTrials(s), ...
        nHits(s) * 100.0 / nTrials(s)), ...
        'units', 'normalized', 'VerticalAlignment', 'top', 'fontSize', 9);
  end
  yl1 = ylim(ax(1));
  yl2 = ylim(ax(2));
  newYL = [min(yl1(1), yl2(1)) max(yl1(2), yl2(2))];
  for a = 1:2
    ylim(ax(a), newYL);
    axes(ax(a));
    % set(gca, 'YLimMode', 'manual');  % freeze the y axes
    patch([preStepVF m m preStepVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
    patch([intStartVF intStopVF intStopVF intStartVF], [newYL(1) newYL(1) newYL(2) newYL(2)], ...
      [0.5 0.5 0.5], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); 
    legend([prefLine(a), probeLine(a)], {"0° ±SEM", sprintf("±%d°", probeDirDeg)}, 'location', 'southwest');
  end
end

function [hPatch, hLine] = plotWithConstSEM(x, y, sem, faceColor)
  
  holdState = ishold;
  hold on;
  yUpper = y + sem;
  yLower = y - sem;
  xPatch = [x(:)', fliplr(x(:)')];
  yPatch = [yUpper(:)', fliplr(yLower(:)')];
  faceAlpha = 0.15;
  hPatch = patch(xPatch, yPatch, [0.7 0.7 0.7], 'faceColor', faceColor, 'FaceAlpha', faceAlpha, ...
    'EdgeColor', 'none');
  hLine = plot(x, y, '-' , 'color', faceColor, 'lineWidth', 1.5);
  uistack(hLine, 'top');
 
  if ~holdState
    hold off;
  end
end

