function T = compareBetaAndKernelScales()
% compareBetaAndKernelScales
% Compare pooled INC beta scales with pooled kernel scales by probe offset.

baseFolder = domainFolder(mfilename('fullpath'));

kernelPath = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', ...
  'IDR_acrossOffsetSummary.mat');

betaPath = fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries', ...
  'IDR_acrossOffsetBetaSummary.mat');

K = load(kernelPath, 'acrossOffsetSummary');
B = load(betaPath, 'betaSummary');

kernelSummary = K.acrossOffsetSummary;
betaSummary = B.betaSummary;

kernelOffsets = [kernelSummary.empirical.probeOffsetDeg]';
kernelScale = [kernelSummary.empirical.pooledScale]';

kernelVar = kernelSummary.bootstrap.fitBootstrap.offsetFitVar(:);
kernelSE = sqrt(kernelVar);

betaOffsets = [betaSummary.offsetFits.probeOffsetDeg]';
betaScale = [betaSummary.offsetFits.scale]';
betaSE = [betaSummary.offsetFits.scaleHessianSE]';

commonOffsets = intersect(kernelOffsets, betaOffsets);

n = numel(commonOffsets);
T = table( ...
  commonOffsets, ...
  nan(n,1), nan(n,1), ...
  nan(n,1), nan(n,1), ...
  'VariableNames', { ...
    'probeOffsetDeg', ...
    'kernelScale', 'kernelSE', ...
    'betaScale', 'betaSE'});

for i = 1:n
  d = commonOffsets(i);

  ik = find(kernelOffsets == d, 1);
  ib = find(betaOffsets == d, 1);

  T.kernelScale(i) = kernelScale(ik);
  T.kernelSE(i) = kernelSE(ik);

  T.betaScale(i) = betaScale(ib);
  T.betaSE(i) = betaSE(ib);
end

T.difference = T.betaScale - T.kernelScale;

fprintf('\nBeta versus kernel scales\n');
disp(T);

fprintf('Correlation across offsets: %.4f\n', ...
  corr(T.kernelScale, T.betaScale, 'Rows', 'complete'));

% ---- Scatterplot with uncertainty in both dimensions ----
figure('Color', 'w', 'Position', [100 100 600 560]);
hold on;

for i = 1:height(T)
  errorbar(T.kernelScale(i), T.betaScale(i), ...
    T.betaSE(i), T.betaSE(i), ...
    T.kernelSE(i), T.kernelSE(i), ...
    'o', ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', 'auto', ...
    'MarkerSize', 7, ...
    'CapSize', 0);

  text(T.kernelScale(i), T.betaScale(i), ...
    sprintf('  %g°', T.probeOffsetDeg(i)), ...
    'VerticalAlignment', 'bottom');
end

allValues = [T.kernelScale; T.betaScale];
lo = min([0; allValues]) - 0.05;
hi = max([1; allValues]) + 0.05;

plot([lo hi], [lo hi], 'k--');
xline(0, ':');
yline(0, ':');

xlabel('Kernel pooled scale');
ylabel('Beta pooled scale');
title('Readout scale: beta regression versus kernels');

xlim([lo hi]);
ylim([lo hi]);
axis square;
box off;

% ---- Offset-space comparison ----
figure('Color', 'w', 'Position', [100 100 750 430]);
hold on;

errorbar(T.probeOffsetDeg, T.kernelScale, T.kernelSE, ...
  'o-', 'LineWidth', 1.2, 'MarkerSize', 6);

errorbar(T.probeOffsetDeg, T.betaScale, T.betaSE, ...
  's-', 'LineWidth', 1.2, 'MarkerSize', 6);

yline(0, ':');
xlabel('Probe direction offset (deg)');
ylabel('Normalized readout scale');
title('Kernel and beta readout measurements');
legend({'Kernel scale', 'Beta scale'}, 'Location', 'best');
box off;

end