function plotBetasByOffset(T)
% plotBetasByOffset  Show sessionwise INC preferred and probe betas by offset.

valid = isfinite(T.probeOffsetDeg) & ...
  isfinite(T.betaPref) & isfinite(T.betaProbe);

T = T(valid, :);
offsets = unique(T.probeOffsetDeg);
nOffsets = numel(offsets);

figure('Color', 'w', 'Position', [100 100 1050 500]);
hold on;

jitterWidth = 1.5;

hPref = gobjects(1);
hProbe = gobjects(1);

for k = 1:nOffsets
  offset = offsets(k);
  idx = T.probeOffsetDeg == offset;
  n = sum(idx);

  % Deterministic horizontal spreading within each offset.
  if n == 1
    jitter = 0;
  else
    jitter = linspace(-jitterWidth, jitterWidth, n)';
  end

  xPref  = offset - 1.8 + jitter;
  xProbe = offset + 1.8 + jitter;

  hp = errorbar(xPref, T.betaPref(idx), T.betaPrefSE(idx), ...
    'o', 'LineStyle', 'none', 'MarkerSize', 4, ...
    'CapSize', 0);

  hq = errorbar(xProbe, T.betaProbe(idx), T.betaProbeSE(idx), ...
    'o', 'LineStyle', 'none', 'MarkerSize', 4, ...
    'CapSize', 0);

  if k == 1
    hPref = hp;
    hProbe = hq;
  end

  % Offset summaries.
  prefMean = mean(T.betaPref(idx), 'omitnan');
  probeMean = mean(T.betaProbe(idx), 'omitnan');

  plot(offset - 1.8, prefMean, 's', ...
    'MarkerFaceColor', 'auto', 'MarkerSize', 8);
  plot(offset + 1.8, probeMean, 's', ...
    'MarkerFaceColor', 'auto', 'MarkerSize', 8);
end

yline(0, ':');
xlabel('Probe direction offset (deg)');
ylabel('Regression beta per % coherence');
title('INC preferred and probe betas by offset');

xticks(offsets);
xlim([min(offsets) - 8, max(offsets) + 8]);

legend([hPref hProbe], ...
  {'Preferred beta', 'Probe beta'}, ...
  'Location', 'best');

box off;
end