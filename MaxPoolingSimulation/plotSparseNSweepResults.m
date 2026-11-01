function figures = plotSparseNSweepResults(results)
% PLOTSPARSENSWEEPRESULTS  Main gain, operating-point, map, and gate figures.

cfg = results.config;
summary = results.summaryTable;
operating = results.operatingTable;
plotNs = intersect(cfg.PlotNValues, cfg.NValues, 'stable');
colors = turbo(max(numel(plotNs),2));

figures = gobjects(4,1);

figures(1) = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Sparse max-pooling gain curves');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
panels = {{'Fixed','Gain','Absolute gain: fixed step'}, ...
          {'Fixed','NormalizedGain','Normalized gain: fixed step'}, ...
          {'Matched','Gain','Absolute gain: 75% matched'}, ...
          {'Matched','NormalizedGain','Normalized gain: 75% matched'}};
for p = 1:4
    nexttile; hold on;
    for k = 1:numel(plotNs)
        take = summary.N == plotNs(k) & summary.Mode == panels{p}{1};
        rows = summary(take,:);
        [~,order] = sort(rows.OffsetDeg);
        rows = rows(order,:);
        plot(rows.OffsetDeg, rows.(panels{p}{2}), '-o', ...
            'Color', colors(k,:), 'LineWidth', 1.1, 'MarkerSize', 3, ...
            'DisplayName', sprintf('n=%d',plotNs(k)));
    end
    yline(0,':','Color',[0.4 0.4 0.4]);
    if strcmp(panels{p}{2},'NormalizedGain')
        yline(1,'--','Color',[0.5 0.5 0.5]);
        tuningOffset = linspace(0,180,361);
        tuningProfile = centeredMTSensitivity(tuningOffset, ...
            cfg.TuningSigmaDeg,cfg.CenterTuning);
        tuningAtDrift = centeredMTSensitivity(0,cfg.TuningSigmaDeg,cfg.CenterTuning);
        plot(tuningOffset,tuningProfile/tuningAtDrift,'k--','LineWidth',1.3, ...
            'DisplayName','Centered MT sensitivity');
    end
    xlim([0 180]); xlabel('Circular offset from drift (deg)');
    ylabel(strrep(panels{p}{2},'NormalizedGain','Gain / drift gain'));
    title(panels{p}{3}); grid on; box off;
end
legend('Location','eastoutside');
sgtitle(sprintf('Sparse aligned mechanisms; sigma_T = %.1f deg',cfg.TuningSigmaDeg));

figures(2) = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Sparse max-pooling operating points');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for mode = ["Fixed","Matched"]
    rows = operating(operating.Mode==mode,:);
    plot(rows.N,rows.Performance,'-o','LineWidth',1.5,'DisplayName',mode);
end
yline(cfg.TargetPerformance,'k:','Target');
xlabel('Number of mechanisms/streams, n'); ylabel('Proportion correct');
ylim([0.45 1]); grid on; box off; legend('Location','best');
title('Behavioral operating point');
nexttile; hold on;
for mode = ["Fixed","Matched"]
    rows = operating(operating.Mode==mode,:);
    plot(rows.N,rows.StepCoherence,'-o','LineWidth',1.5,'DisplayName',mode);
end
xlabel('Number of mechanisms/streams, n'); ylabel('Deterministic drift step');
grid on; box off; legend('Location','best'); title('Step required by operating mode');

figures(3) = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Sparse max-pooling gain maps');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for p = 1:2
    nexttile;
    modeName = ["Fixed","Matched"]; modeName = modeName(p);
    rows = summary(summary.Mode==modeName,:);
    scatter(rows.OffsetDeg,rows.N,55,rows.NormalizedGain,'filled');
    xlabel('Circular offset from drift (deg)'); ylabel('n');
    xlim([0 180]); colorbar; colormap(turbo); box off; grid on;
    title(sprintf('%s: normalized gain',modeName));
end

figures(4) = figure('Color','w','Visible',cfg.Visible, ...
    'Name','Close-to-max diagnostics');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
d = results.diagnostic;
valid = d.binCounts > 25 & isfinite(d.meanCounterfactualInfluence);
plot(d.gapBinCenters(valid),d.meanCounterfactualInfluence(valid),'-o', ...
    'LineWidth',1.5,'MarkerSize',4);
xline(0,'k:','Winner boundary');
xlabel('Competitive gap: R_m - maximum of other mechanisms', ...
    'Interpreter','none');
ylabel('Counterfactual change in max / perturbation');
ylim([-0.02 1.02]); grid on; box off;
title(sprintf('Mechanism gate (n=%d)',d.n));
nexttile;
rows = summary(summary.Mode==string(cfg.DiagnosticMode) & summary.OffsetDeg>0,:);
scatter(rows.AlignedMechanismWinProbability,rows.NormalizedGain,55, ...
    rows.OffsetDeg,'filled');
xlabel('Aligned mechanism win probability'); ylabel('Normalized stream gain');
colorbar; colormap(turbo); grid on; box off;
title('Recovered gain and winner access');
end
