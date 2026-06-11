function fig = plotKernelVsScalarScale(probeSummary, kernelScaleTable)

% probeSummary needs:
%   probeTag or probeDirDeg
%   ratioOfMeans_changeProbeOverPref_raw_boot
%   ratioOfMeans_changeProbeOverPref_raw_ciLow
%   ratioOfMeans_changeProbeOverPref_raw_ciHigh
%
% kernelScaleTable needs:
%   probeDirDeg
%   scale
% optionally:
%   ciLow
%   ciHigh

if ismember('probeDirDeg', probeSummary.Properties.VariableNames)
  xScalar = probeSummary.probeDirDeg;
else
  xScalar = nan(height(probeSummary), 1);
  for i = 1:height(probeSummary)
    xScalar(i) = sscanf(char(probeSummary.probeTag(i)), 'Probe%d');
  end
end

yScalar = probeSummary.ratioOfMeans_changeProbeOverPref_raw_boot;
loScalar = probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciLow;
hiScalar = probeSummary.ratioOfMeans_changeProbeOverPref_raw_ciHigh;

[xScalar, order] = sort(xScalar);
yScalar = yScalar(order);
loScalar = loScalar(order);
hiScalar = hiScalar(order);

errLow = yScalar - loScalar;
errHigh = hiScalar - yScalar;

fig = figure('Color', 'w', 'Position', [100 100 750 450]);
hold on;

errorbar(xScalar, yScalar, errLow, errHigh, '-o', ...
  'LineWidth', 1.5, ...
  'MarkerFaceColor', 'auto', ...
  'DisplayName', 'Scalar regression');

if all(ismember({'probeDirDeg', 'scale'}, kernelScaleTable.Properties.VariableNames))

  [xKernel, orderK] = sort(kernelScaleTable.probeDirDeg);
  yKernel = kernelScaleTable.scale(orderK);

  if all(ismember({'ciLow', 'ciHigh'}, kernelScaleTable.Properties.VariableNames))
    loK = kernelScaleTable.ciLow(orderK);
    hiK = kernelScaleTable.ciHigh(orderK);
    errorbar(xKernel, yKernel, yKernel - loK, hiK - yKernel, '-s', ...
      'LineWidth', 1.5, ...
      'MarkerFaceColor', 'auto', ...
      'DisplayName', 'Kernel scale');
  else
    plot(xKernel, yKernel, '-s', ...
      'LineWidth', 1.5, ...
      'MarkerFaceColor', 'auto', ...
      'DisplayName', 'Kernel scale');
  end
end

yline(0, ':');
yline(1, ':');

xlabel('Probe offset (deg)');
ylabel('Probe / preferred scale');
title('Kernel scale versus scalar-regression scale');
legend('Location', 'best');
box off;

end