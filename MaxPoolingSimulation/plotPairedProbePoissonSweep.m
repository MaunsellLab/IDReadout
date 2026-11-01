function handles = plotPairedProbePoissonSweep(results)
% PLOTPAIREDPROBEPOISSONSWEEP  Base-drift, gain, threshold, and winner checks.

T=results.conditionTable; cfg=results.config;
bases=unique(T.BaseDriftCoherence,'stable');
modes=["Independent","Yoked"];
colors=lines(numel(bases));
handles=gobjects(4,1);

handles(1)=figure('Visible',cfg.Visible,'Color','w','Name','Poisson gains');
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    ax=nexttile; hold(ax,'on');
    for bb=1:numel(bases)
        rows=T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb);
        S=sortrows(T(rows,:), 'OffsetDeg');
        plot(ax,S.OffsetDeg,S.DriftGain,'-o','Color',colors(bb,:), ...
            'DisplayName',sprintf('Drift, base %.3g%%',bases(bb)));
        plot(ax,S.OffsetDeg,S.ProbeGainPerStream,'--s','Color',colors(bb,:), ...
            'DisplayName',sprintf('Probe, base %.3g%%',bases(bb)));
    end
    H=results.homoscedasticReferenceTable;
    if ~isempty(H)
        rows=H.Yoking==modes(mm);
        S=sortrows(H(rows,:), 'OffsetDeg');
        plot(ax,S.OffsetDeg,S.DriftGain,':','Color',[0.25 0.25 0.25], ...
            'LineWidth',1.5,'DisplayName','Drift, homoscedastic');
        plot(ax,S.OffsetDeg,S.ProbeGainPerStream,'-.','Color',[0.25 0.25 0.25], ...
            'LineWidth',1.5,'DisplayName','Probe, homoscedastic');
    end
    yline(ax,0,':k'); xlabel(ax,'Paired-probe offset (deg)');
    ylabel(ax,'Recovered probit gain'); title(ax,modes(mm)); grid(ax,'on');
    if mm==2, legend(ax,'Location','best'); end
end
title(tl,'Response-dependent pooled-Poisson variance');

handles(2)=figure('Visible',cfg.Visible,'Color','w','Name','Poisson operating points');
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    ax=nexttile; hold(ax,'on');
    for bb=1:numel(bases)
        rows=T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb);
        S=sortrows(T(rows,:), 'OffsetDeg');
        plot(ax,S.OffsetDeg,100*S.Performance,'-o','Color',colors(bb,:), ...
            'DisplayName',sprintf('20%% performance, base %.3g%%',bases(bb)));
        plot(ax,S.OffsetDeg,S.Threshold75,'--s','Color',colors(bb,:), ...
            'DisplayName',sprintf('75%% threshold, base %.3g%%',bases(bb)));
    end
    yline(ax,75,':k','75% correct'); yline(ax,20,':','20% increment');
    xlabel(ax,'Paired-probe offset (deg)'); ylabel(ax,'Percent or coherence %');
    title(ax,modes(mm)); grid(ax,'on');
    if mm==2, legend(ax,'Location','best'); end
end
title(tl,'Fixed-increment performance and independently measured threshold');

handles(3)=figure('Visible',cfg.Visible,'Color','w','Name','Poisson winners');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    for bb=1:numel(bases)
        ax=nexttile; hold(ax,'on');
        rows=T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb);
        S=sortrows(T(rows,:), 'OffsetDeg');
        plot(ax,S.OffsetDeg,S.ChangedDriftWinner,'-o','DisplayName','Changed drift');
        plot(ax,S.OffsetDeg,S.ChangedProbeWinner,'-s','DisplayName','Changed probe');
        plot(ax,S.OffsetDeg,S.ChangedOtherWinner,'-^','DisplayName','Changed other');
        plot(ax,S.OffsetDeg,S.NoChangeOtherWinner,'--','LineWidth',1.5, ...
            'DisplayName','No-change other');
        ylim(ax,[0 1]); xlabel(ax,'Offset (deg)'); ylabel(ax,'Winner fraction');
        title(ax,sprintf('%s, base %.3g%%',modes(mm),bases(bb))); grid(ax,'on');
        if mm==1 && bb==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Location of the hard-max winner');

handles(4)=figure('Visible',cfg.Visible,'Color','w','Name','Poisson rate diagnostics');
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile; hold(ax,'on');
for bb=1:numel(bases)
    rows=T.Yoking=="Independent" & T.BaseDriftCoherence==bases(bb);
    S=sortrows(T(rows,:), 'OffsetDeg');
    plot(ax,S.OffsetDeg,S.PreRateDrift,'-o','Color',colors(bb,:), ...
        'DisplayName',sprintf('Drift, base %.3g%%',bases(bb)));
    plot(ax,S.OffsetDeg,S.PreRateProbe,'--s','Color',colors(bb,:), ...
        'DisplayName',sprintf('Probe, base %.3g%%',bases(bb)));
end
xlabel(ax,'Offset (deg)'); ylabel(ax,'Prestep rate (sp/s)'); grid(ax,'on');
legend(ax,'Location','best'); title(ax,'Absolute rates that set count variance');
ax=nexttile;
bar(ax,bases,arrayfun(@(b)max(T.ChangedClampFraction( ...
    T.BaseDriftCoherence==b)),bases));
xlabel(ax,'Base drift coherence (%)'); ylabel(ax,'Maximum clamp fraction');
title(ax,sprintf('Calibrated N_{eff} = %.3g', ...
    results.poolCalibration.effectivePoolSize)); grid(ax,'on');
title(tl,'Rate and positivity diagnostics');
end
