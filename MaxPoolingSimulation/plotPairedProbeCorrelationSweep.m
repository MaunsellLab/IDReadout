function figures = plotPairedProbeCorrelationSweep(results)
% PLOTPAIREDPROBECORRELATIONSWEEP  Summary figures across noise smoothness.

cfg=results.config;
T=results.conditionTable;
C=results.noiseCalibrationTable;
lengths=cfg.CorrelationLengthsDeg(:)';
modes=["Independent","Yoked"];
colors=turbo(numel(lengths));
figures=gobjects(5,1);

figures(1)=figure('Color','w','Visible',cfg.Visible, ...
    'Name','Correlation sweep gains');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for mm=1:2
    nexttile; hold on;
    for ll=1:numel(lengths)
        rows=T(T.Yoking==modes(mm) & ...
            T.InternalCorrelationLengthDeg==lengths(ll),:);
        rows=sortrows(rows,'OffsetDeg');
        plot(rows.OffsetDeg,rows.NormalizedProbePerStream,'-o', ...
            'Color',colors(ll,:),'LineWidth',1.3,'MarkerSize',3, ...
            'DisplayName',sprintf('L=%g deg, rank=%.1f',lengths(ll), ...
            rows.EffectiveNoiseRank(1)));
    end
    yline(0,'k:'); xlim([0 180]); grid on; box off;
    xlabel('Paired-probe offset (deg)'); ylabel('Probe / drift gain');
    title(modes(mm));
end
legend('Location','eastoutside');
sgtitle('Effect of internal-noise correlation length on recovered gain');

figures(2)=figure('Color','w','Visible',cfg.Visible, ...
    'Name','Correlation sweep operating points');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for mm=1:2
    nexttile(mm); hold on;
    for ll=1:numel(lengths)
        rows=T(T.Yoking==modes(mm) & ...
            T.InternalCorrelationLengthDeg==lengths(ll),:);
        rows=sortrows(rows,'OffsetDeg');
        plot(rows.OffsetDeg,rows.Threshold75,'-','Color',colors(ll,:), ...
            'LineWidth',1.2,'DisplayName',sprintf('L=%g',lengths(ll)));
    end
    yline(cfg.PairedConfig.FixedStepCoherence,'k:'); grid on; box off;
    xlabel('Offset (deg)'); ylabel('75% threshold'); title(modes(mm));

    nexttile(mm+2); hold on;
    for ll=1:numel(lengths)
        rows=T(T.Yoking==modes(mm) & ...
            T.InternalCorrelationLengthDeg==lengths(ll),:);
        rows=sortrows(rows,'OffsetDeg');
        plot(rows.OffsetDeg,rows.Performance,'-','Color',colors(ll,:), ...
            'LineWidth',1.2);
    end
    yline(cfg.PairedConfig.TargetPerformance,'k:'); grid on; box off;
    xlabel('Offset (deg)'); ylabel(sprintf('Performance at %.1f%%', ...
        cfg.PairedConfig.FixedStepCoherence));
end
sgtitle('Fixed-noise operating points');

[~,order]=sort(C.correlationLengthDeg);
C=C(order,:);
figures(3)=figure('Color','w','Visible',cfg.Visible, ...
    'Name','Correlation sweep noise calibration');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile;
plot(C.correlationLengthDeg,C.internalNoiseSD,'-o','LineWidth',1.5);
grid on; box off; xlabel('Correlation length (deg)');
ylabel('Calibrated internal-noise SD');
title(sprintf('Noise required for %.0f%% at %.1f%%', ...
    100*cfg.PairedConfig.TargetPerformance,cfg.PairedConfig.FixedStepCoherence));
nexttile;
plot(C.correlationLengthDeg,C.effectiveNoiseRank,'-o','LineWidth',1.5);
grid on; box off; xlabel('Correlation length (deg)');
ylabel('Effective covariance rank');
title('After circular-mean removal');

selectedOffsets=intersect([45 90 120 180],unique(T.OffsetDeg)','stable');
figures(4)=figure('Color','w','Visible',cfg.Visible, ...
    'Name','Correlation sweep selected offsets');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
offsetColors=lines(max(numel(selectedOffsets),1));
for mm=1:2
    nexttile; hold on;
    for oo=1:numel(selectedOffsets)
        rows=T(T.Yoking==modes(mm) & T.OffsetDeg==selectedOffsets(oo),:);
        rows=sortrows(rows,'EffectiveNoiseRank');
        plot(rows.EffectiveNoiseRank,rows.NormalizedProbePerStream,'-o', ...
            'Color',offsetColors(oo,:),'LineWidth',1.4, ...
            'DisplayName',sprintf('%g deg',selectedOffsets(oo)));
    end
    yline(0,'k:'); grid on; box off;
    xlabel('Effective covariance rank'); ylabel('Probe / drift gain');
    title(modes(mm));
end
legend('Location','best');
sgtitle('Selected gain values versus effective noise dimensionality');

winnerOffsets=intersect([45 120 180],unique(T.OffsetDeg)','stable');
figures(5)=figure('Color','w','Visible',cfg.Visible, ...
    'Name','Correlation sweep other winners');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for mm=1:2
    for patch=1:2
        nexttile((patch-1)*2+mm); hold on;
        if patch==1, variable='ChangedOtherWinner'; patchName='Changed';
        else, variable='NoChangeOtherWinner'; patchName='No-change'; end
        for oo=1:numel(winnerOffsets)
            rows=T(T.Yoking==modes(mm) & T.OffsetDeg==winnerOffsets(oo),:);
            rows=sortrows(rows,'EffectiveNoiseRank');
            plot(rows.EffectiveNoiseRank,rows.(variable),'-o', ...
                'Color',offsetColors(oo,:),'LineWidth',1.4, ...
                'DisplayName',sprintf('%g deg',winnerOffsets(oo)));
        end
        ylim([0 1]); grid on; box off;
        xlabel('Effective covariance rank'); ylabel('Other-winner fraction');
        title(sprintf('%s: %s patch',modes(mm),patchName));
    end
end
legend('Location','best');
sgtitle('Categorical other-winner diagnostic');
end
