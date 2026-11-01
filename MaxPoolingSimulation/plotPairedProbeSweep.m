function figures = plotPairedProbeSweep(results)
% PLOTPAIREDPROBESWEEP  Gains, winners, operating points, and decompositions.

cfg=results.config; T=results.conditionTable; D=results.decompositionTable;
modes=["Independent","Yoked"]; colors=lines(2);
figures=gobjects(4,1);

figures(1)=figure('Color','w','Visible',cfg.Visible,'Name','Paired probe gains');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
variables={'DriftGain','ProbeGainPerStream','NormalizedProbePerStream','ProbeGainTotal'};
titles={'Drift gain','Probe gain: per-stream equivalent', ...
    'Probe / drift gain','Total paired-probe gain'};
for p=1:4
    nexttile; hold on;
    for m=1:2
        rows=T(T.Yoking==modes(m),:);
        plot(rows.OffsetDeg,rows.(variables{p}),'-o','Color',colors(m,:), ...
            'LineWidth',1.5,'MarkerSize',4,'DisplayName',modes(m));
    end
    yline(0,'k:'); xlim([0 180]); grid on; box off;
    xlabel('Paired-probe offset (deg)'); ylabel(titles{p}); title(titles{p});
end
legend('Location','best');
sgtitle(sprintf('Delta-coherence paired probes: fixed step %.1f, sigma_T %.1f deg', ...
    cfg.FixedStepCoherence,cfg.TuningSigmaDeg));

figures(2)=figure('Color','w','Visible',cfg.Visible,'Name','Paired probe winners');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for m=1:2
    nexttile; hold on; rows=T(T.Yoking==modes(m),:);
    plot(rows.OffsetDeg,rows.ChangedDriftWinner,'-','LineWidth',1.5,'DisplayName','Drift-associated');
    plot(rows.OffsetDeg,rows.ChangedProbeWinner,'-','LineWidth',1.5,'DisplayName','Probe-associated');
    plot(rows.OffsetDeg,rows.ChangedOtherWinner,'-','LineWidth',1.5,'DisplayName','Other/intermediate');
    xlabel('Paired-probe offset (deg)'); ylabel('Changed-patch winner fraction');
    ylim([0 1]); grid on; box off; title(modes(m)); legend('Location','best');
end

figures(3)=figure('Color','w','Visible',cfg.Visible,'Name','Paired probe operating points');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for m=1:2
    rows=T(T.Yoking==modes(m),:);
    plot(rows.OffsetDeg,rows.Threshold75,'-o','Color',colors(m,:), ...
        'LineWidth',1.5,'DisplayName',modes(m));
end
yline(cfg.FixedStepCoherence,'k:',sprintf('%.1f%% reference',cfg.FixedStepCoherence));
xlabel('Paired-probe offset (deg)'); ylabel('75% drift-increment threshold'); grid on; box off; legend('Location','best');
title('Threshold with fixed noise');
nexttile; hold on;
for m=1:2
    rows=T(T.Yoking==modes(m),:);
    plot(rows.OffsetDeg,rows.Performance,'-o','Color',colors(m,:), ...
        'LineWidth',1.5,'DisplayName',modes(m));
end
yline(cfg.TargetPerformance,'k:'); ylim([0.45 1]);
xlabel('Paired-probe offset (deg)'); ylabel('Proportion correct'); grid on; box off;
title(sprintf('Performance at fixed %.1f%% increment',cfg.FixedStepCoherence));

figures(4)=figure('Color','w','Visible',cfg.Visible,'Name','Stream influence decomposition');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for m=1:2
    for s=1:2
        nexttile((s-1)*2+m); hold on;
        stream=["ProbePerStream","Drift"]; stream=stream(s);
        rows=D(D.Yoking==modes(m) & D.Stream==stream,:);
        plot(rows.OffsetDeg,rows.DecisionFixedWinner,'-','LineWidth',1.5,'DisplayName','Fixed winner');
        plot(rows.OffsetDeg,rows.DecisionWinnerSwitch,'-','LineWidth',1.5,'DisplayName','Winner switch');
        plot(rows.OffsetDeg,rows.DecisionTotal,'k--','LineWidth',1.2,'DisplayName','Total evidence derivative');
        yline(0,'k:'); xlim([0 180]); grid on; box off;
        xlabel('Paired-probe offset (deg)'); ylabel('Evidence derivative');
        title(sprintf('%s: %s',modes(m),stream));
        if m==2 && s==1, legend('Location','best'); end
    end
end
end
