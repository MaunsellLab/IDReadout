function handles = plotPairedProbeTemporalMaxSweep(results)
% PLOTPAIREDPROBETEMPORALMAXSWEEP  Gains, thresholds, transitions, mechanisms.

T=results.conditionTable; D=results.decompositionTable;
W=results.winnerDirectionTable; cfg=results.config;
bases=unique(T.BaseDriftCoherence,'stable'); modes=["Independent","Yoked"];
colors=lines(numel(bases)); handles=gobjects(6,1);

handles(1)=figure('Visible',cfg.Visible,'Color','w','Name','Temporal-max gains');
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    ax=nexttile; hold(ax,'on');
    for bb=1:numel(bases)
        S=sortrows(T(T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb),:), ...
            'OffsetDeg');
        plot(ax,S.OffsetDeg,S.DriftGain,'-o','Color',colors(bb,:), ...
            'DisplayName',sprintf('Drift, base %.3g%%',bases(bb)));
        plot(ax,S.OffsetDeg,S.ProbeGainPerStream,'--s','Color',colors(bb,:), ...
            'DisplayName',sprintf('Probe, base %.3g%%',bases(bb)));
    end
    yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
    xlabel(ax,'Paired-probe offset (deg)'); ylabel(ax,'Recovered probit gain');
    title(ax,modes(mm)); if mm==2, legend(ax,'Location','best'); end
end
title(tl,'Separate prestep and step maxima: recovered gains');

handles(2)=figure('Visible',cfg.Visible,'Color','w','Name','Temporal-max operating points');
tl=tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    ax=nexttile; hold(ax,'on');
    for bb=1:numel(bases)
        S=sortrows(T(T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb),:), ...
            'OffsetDeg');
        plot(ax,S.OffsetDeg,100*S.Performance,'-o','Color',colors(bb,:), ...
            'DisplayName',sprintf('20%% performance, base %.3g%%',bases(bb)));
        plot(ax,S.OffsetDeg,S.Threshold75,'--s','Color',colors(bb,:), ...
            'DisplayName',sprintf('75%% threshold, base %.3g%%',bases(bb)));
    end
    yline(ax,75,':k','HandleVisibility','off');
    yline(ax,20,':','HandleVisibility','off'); grid(ax,'on');
    xlabel(ax,'Paired-probe offset (deg)'); ylabel(ax,'Percent or coherence %');
    title(ax,modes(mm)); if mm==2, legend(ax,'Location','best'); end
end
title(tl,'Fixed-increment performance and independently measured threshold');

handles(3)=figure('Visible',cfg.Visible,'Color','w','Name','Temporal winner dynamics');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    for bb=1:numel(bases)
        ax=nexttile; hold(ax,'on');
        S=sortrows(T(T.Yoking==modes(mm) & T.BaseDriftCoherence==bases(bb),:), ...
            'OffsetDeg');
        plot(ax,S.OffsetDeg,S.ChangedTemporalSameWinner,'-o', ...
            'DisplayName','Changed: identical mechanism');
        plot(ax,S.OffsetDeg,S.NoChangeTemporalSameWinner,'--s', ...
            'DisplayName','No-change: identical mechanism');
        plot(ax,S.OffsetDeg,S.ChangedDriftWinner,'-','DisplayName','Changed step: drift');
        plot(ax,S.OffsetDeg,S.ChangedProbeWinner,'-','DisplayName','Changed step: probe');
        ylim(ax,[0 1]); grid(ax,'on'); xlabel(ax,'Offset (deg)');
        ylabel(ax,'Trial fraction');
        title(ax,sprintf('%s, base %.3g%%',modes(mm),bases(bb)));
        if mm==1 && bb==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Prestep-to-step winner stability and step-winner identity');

handles(4)=figure('Visible',cfg.Visible,'Color','w','Name','Temporal gain decomposition');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
streams=["Drift","ProbePerStream"];
for mm=1:2
    for ss=1:2
        ax=nexttile((ss-1)*2+mm); hold(ax,'on');
        for bb=1:numel(bases)
            S=sortrows(D(D.Yoking==modes(mm) & D.Stream==streams(ss) & ...
                D.BaseDriftCoherence==bases(bb),:),'OffsetDeg');
            plot(ax,S.OffsetDeg,S.TemporalSameWinner,'-','Color',colors(bb,:), ...
                'DisplayName',sprintf('Same, base %.3g%%',bases(bb)));
            plot(ax,S.OffsetDeg,S.TemporalDifferentWinner,'--','Color',colors(bb,:), ...
                'DisplayName',sprintf('Different, base %.3g%%',bases(bb)));
            plot(ax,S.OffsetDeg,S.PerturbationWinnerSwitch,':','Color',colors(bb,:), ...
                'DisplayName',sprintf('Perturb switch, base %.3g%%',bases(bb)));
        end
        yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
        xlabel(ax,'Offset (deg)'); ylabel(ax,'Evidence derivative');
        title(ax,sprintf('%s: %s',modes(mm),streams(ss)));
        if mm==2 && ss==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Stream influence by temporal winner relationship');

handles(5)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Winner-direction distributions');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases)
    for mm=1:2
        ax=nexttile((bb-1)*2+mm); hold(ax,'on');
        S=selectWinnerDirectionCondition(W,bases(bb),modes(mm), ...
            cfg.WinnerDirectionPlotOffsetDeg);
        stairs(ax,S.WinnerDirectionDeg,S.ChangedWinnerFraction,'-', ...
            'DisplayName','Changed patch');
        stairs(ax,S.WinnerDirectionDeg,S.NoChangeWinnerFraction,'--', ...
            'DisplayName','No-change patch');
        stairs(ax,S.WinnerDirectionDeg,S.PatchAverageWinnerFraction,':', ...
            'LineWidth',1.5,'DisplayName','Patch average');
        xlim(ax,[-180 180]); grid(ax,'on');
        xlabel(ax,'Step-winner direction relative to drift (deg)');
        ylabel(ax,sprintf('Probability per %.3g-deg bin', ...
            cfg.WinnerDirectionBinWidthDeg));
        title(ax,sprintf('%s, base %.3g%%, offset %.3g deg', ...
            modes(mm),bases(bb),S.OffsetDeg(1)));
        if bb==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Distribution of step winners across the dense direction bank');

handles(6)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Winner-direction signed contributions');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases)
    for mm=1:2
        ax=nexttile((bb-1)*2+mm); hold(ax,'on');
        S=selectWinnerDirectionCondition(W,bases(bb),modes(mm), ...
            cfg.WinnerDirectionPlotOffsetDeg);
        stairs(ax,S.WinnerDirectionDeg,S.DriftWinnerWeightedContribution,'-', ...
            'LineWidth',1.4,'DisplayName','Drift');
        stairs(ax,S.WinnerDirectionDeg, ...
            S.ProbeWinnerWeightedContributionPerStream,'--', ...
            'LineWidth',1.4,'DisplayName','Probe per stream');
        yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[-180 180]);
        grid(ax,'on');
        xlabel(ax,'Step-winner direction relative to drift (deg)');
        ylabel(ax,'Winner probability x local sensitivity');
        title(ax,sprintf('%s, base %.3g%%: sum D=%.3f, P=%.3f', ...
            modes(mm),bases(bb),S.DriftFixedWinnerTotal(1), ...
            S.ProbeFixedWinnerTotalPerStream(1)));
        if bb==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,sprintf(['Signed fixed-winner contribution by direction ' ...
    '(offset nearest %.3g deg)'],cfg.WinnerDirectionPlotOffsetDeg));
end

function S=selectWinnerDirectionCondition(W,base,mode,targetOffset)
available=unique(W.OffsetDeg(W.BaseDriftCoherence==base & W.Yoking==mode));
[~,which]=min(abs(available-targetOffset)); offset=available(which);
S=sortrows(W(W.BaseDriftCoherence==base & W.Yoking==mode & ...
    W.OffsetDeg==offset,:),'WinnerDirectionDeg');
end
