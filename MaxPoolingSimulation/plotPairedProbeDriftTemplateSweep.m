function handles = plotPairedProbeDriftTemplateSweep(results)
% PLOTPAIREDPROBEDRIFTTEMPLATESWEEP  Alpha, gains, winners, and contributions.

T=results.conditionTable; W=results.winnerDirectionTable; cfg=results.config;
bases=unique(T.BaseDriftCoherence,'stable');
offsets=unique(T.OffsetDeg,'stable'); modes=["Independent","Yoked"];
colors=lines(numel(bases)); offsetColors=lines(numel(offsets));
handles=gobjects(6,1);

handles(1)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template gains');
tl=tiledlayout(numel(offsets),2,'Padding','compact','TileSpacing','compact');
for oo=1:numel(offsets)
    for mm=1:2
        ax=nexttile((oo-1)*2+mm); hold(ax,'on');
        for bb=1:numel(bases)
            S=sortrows(T(T.OffsetDeg==offsets(oo) & T.Yoking==modes(mm) & ...
                T.BaseDriftCoherence==bases(bb),:),'CancellationFraction');
            plot(ax,S.CancellationFraction,S.DriftGain,'-o', ...
                'Color',colors(bb,:),'DisplayName',sprintf('Drift, base %.3g%%',bases(bb)));
            plot(ax,S.CancellationFraction,S.ProbeGainPerStream,'--s', ...
                'Color',colors(bb,:),'DisplayName',sprintf('Probe, base %.3g%%',bases(bb)));
        end
        yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
        xlabel(ax,'Drift-template cancellation fraction');
        ylabel(ax,'Recovered probit gain');
        title(ax,sprintf('%s, offset %.3g deg',modes(mm),offsets(oo)));
        if oo==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Recovered gains as the perceived drift template is removed');

handles(2)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template operating points');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    for rr=1:2
        ax=nexttile((rr-1)*2+mm); hold(ax,'on');
        for oo=1:numel(offsets)
            for bb=1:numel(bases)
                S=sortrows(T(T.OffsetDeg==offsets(oo) & T.Yoking==modes(mm) & ...
                    T.BaseDriftCoherence==bases(bb),:),'CancellationFraction');
                style='-'; if bb==2, style='--'; end
                if rr==1, value=S.Threshold75; else, value=100*S.Performance; end
                plot(ax,S.CancellationFraction,value,style,'Color',offsetColors(oo,:), ...
                    'DisplayName',sprintf('Offset %.3g, base %.3g%%', ...
                    offsets(oo),bases(bb)));
            end
        end
        if rr==1
            yline(ax,cfg.PairedConfig.FixedStepCoherence,':k', ...
                'HandleVisibility','off'); ylabel(ax,'75% threshold (coherence %)');
        else
            yline(ax,100*cfg.PairedConfig.TargetPerformance,':k', ...
                'HandleVisibility','off'); ylabel(ax,'Performance at 20% (%)');
        end
        grid(ax,'on'); xlabel(ax,'Drift-template cancellation fraction');
        title(ax,modes(mm));
        if rr==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Operating-point stability with one fixed internal-noise SD');

handles(3)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template winner fractions');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases)
    for mm=1:2
        ax=nexttile((bb-1)*2+mm); hold(ax,'on');
        S=selectCondition(T,bases(bb),modes(mm),cfg.WinnerDirectionPlotOffsetDeg);
        plot(ax,S.CancellationFraction,S.ChangedProbeWinner,'-o', ...
            'DisplayName','Changed: probe');
        plot(ax,S.CancellationFraction,S.NoChangeProbeWinner,'--o', ...
            'DisplayName','No-change: probe');
        plot(ax,S.CancellationFraction,S.ChangedDriftWinner,'-s', ...
            'DisplayName','Changed: drift');
        plot(ax,S.CancellationFraction,S.NoChangeDriftWinner,'--s', ...
            'DisplayName','No-change: drift');
        ylim(ax,[0 1]); grid(ax,'on');
        xlabel(ax,'Drift-template cancellation fraction'); ylabel(ax,'Winner fraction');
        title(ax,sprintf('%s, base %.3g%%, offset %.3g deg', ...
            modes(mm),bases(bb),S.OffsetDeg(1)));
        if bb==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Drift- and probe-associated winners reallocated by cancellation');

handles(4)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template fixed-winner totals');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases)
    for mm=1:2
        ax=nexttile((bb-1)*2+mm); hold(ax,'on');
        S=selectWinnerTotals(W,bases(bb),modes(mm), ...
            cfg.WinnerDirectionPlotOffsetDeg);
        plot(ax,S.CancellationFraction,S.DriftFixedWinnerTotal,'-o', ...
            'DisplayName','Drift');
        plot(ax,S.CancellationFraction,S.ProbeFixedWinnerTotalPerStream,'--s', ...
            'DisplayName','Probe per stream');
        yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
        xlabel(ax,'Drift-template cancellation fraction');
        ylabel(ax,'Winner-weighted fixed-winner derivative');
        title(ax,sprintf('%s, base %.3g%%, offset %.3g deg', ...
            modes(mm),bases(bb),S.OffsetDeg(1)));
        if bb==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,'Signed contribution integrated over winner direction');

[plotBase,plotOffset]=nearestBaseOffset(W,cfg.WinnerDirectionPlotBase, ...
    cfg.WinnerDirectionPlotOffsetDeg);
plotFractions=nearestFractions(W,plotBase,plotOffset, ...
    cfg.WinnerDirectionPlotFractions);
handles(5)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template winner-direction distributions');
tl=tiledlayout(numel(plotFractions),2,'Padding','compact','TileSpacing','compact');
for aa=1:numel(plotFractions)
    for mm=1:2
        ax=nexttile((aa-1)*2+mm); hold(ax,'on');
        S=selectDirectionRows(W,plotBase,plotOffset,modes(mm),plotFractions(aa));
        stairs(ax,S.WinnerDirectionDeg,S.ChangedWinnerFraction,'-', ...
            'DisplayName','Changed patch');
        stairs(ax,S.WinnerDirectionDeg,S.NoChangeWinnerFraction,'--', ...
            'DisplayName','No-change patch');
        stairs(ax,S.WinnerDirectionDeg,S.PatchAverageWinnerFraction,':', ...
            'LineWidth',1.4,'DisplayName','Patch average');
        xlim(ax,[-180 180]); grid(ax,'on');
        xlabel(ax,'Step-winner direction relative to drift (deg)');
        ylabel(ax,sprintf('Probability per %.3g-deg bin', ...
            cfg.WinnerDirectionBinWidthDeg));
        title(ax,sprintf('%s, alpha %.3g',modes(mm),plotFractions(aa)));
        if aa==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,sprintf('Winner distributions: base %.3g%%, offset %.3g deg', ...
    plotBase,plotOffset));

handles(6)=figure('Visible',cfg.Visible,'Color','w', ...
    'Name','Drift-template winner-direction contributions');
tl=tiledlayout(numel(plotFractions),2,'Padding','compact','TileSpacing','compact');
for aa=1:numel(plotFractions)
    for mm=1:2
        ax=nexttile((aa-1)*2+mm); hold(ax,'on');
        S=selectDirectionRows(W,plotBase,plotOffset,modes(mm),plotFractions(aa));
        stairs(ax,S.WinnerDirectionDeg,S.DriftWinnerWeightedContribution,'-', ...
            'DisplayName','Drift');
        stairs(ax,S.WinnerDirectionDeg, ...
            S.ProbeWinnerWeightedContributionPerStream,'--', ...
            'DisplayName','Probe per stream');
        yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[-180 180]); grid(ax,'on');
        xlabel(ax,'Step-winner direction relative to drift (deg)');
        ylabel(ax,'Winner probability x local sensitivity');
        title(ax,sprintf('%s, alpha %.3g: sum D=%.3f, P=%.3f', ...
            modes(mm),plotFractions(aa),S.DriftFixedWinnerTotal(1), ...
            S.ProbeFixedWinnerTotalPerStream(1)));
        if aa==1 && mm==1, legend(ax,'Location','best'); end
    end
end
title(tl,sprintf('Signed contributions: base %.3g%%, offset %.3g deg', ...
    plotBase,plotOffset));
end

function S=selectCondition(T,base,mode,targetOffset)
available=unique(T.OffsetDeg(T.BaseDriftCoherence==base & T.Yoking==mode));
[~,which]=min(abs(available-targetOffset)); offset=available(which);
S=sortrows(T(T.BaseDriftCoherence==base & T.Yoking==mode & ...
    T.OffsetDeg==offset,:),'CancellationFraction');
end

function S=selectWinnerTotals(W,base,mode,targetOffset)
R=selectDirectionRowsForOffset(W,base,mode,targetOffset);
[~,first]=unique(R.CancellationFraction,'stable');
S=sortrows(R(first,:),'CancellationFraction');
end

function R=selectDirectionRowsForOffset(W,base,mode,targetOffset)
available=unique(W.OffsetDeg(W.BaseDriftCoherence==base & W.Yoking==mode));
[~,which]=min(abs(available-targetOffset)); offset=available(which);
R=W(W.BaseDriftCoherence==base & W.Yoking==mode & W.OffsetDeg==offset,:);
end

function [base,offset]=nearestBaseOffset(W,targetBase,targetOffset)
bases=unique(W.BaseDriftCoherence); [~,bb]=min(abs(bases-targetBase)); base=bases(bb);
offsets=unique(W.OffsetDeg(W.BaseDriftCoherence==base));
[~,oo]=min(abs(offsets-targetOffset)); offset=offsets(oo);
end

function fractions=nearestFractions(W,base,offset,targets)
available=unique(W.CancellationFraction(W.BaseDriftCoherence==base & ...
    W.OffsetDeg==offset));
fractions=zeros(size(targets));
for k=1:numel(targets)
    [~,which]=min(abs(available-targets(k))); fractions(k)=available(which);
end
fractions=unique(fractions,'stable');
end

function S=selectDirectionRows(W,base,offset,mode,alpha)
S=sortrows(W(W.BaseDriftCoherence==base & W.OffsetDeg==offset & ...
    W.Yoking==mode & W.CancellationFraction==alpha,:),'WinnerDirectionDeg');
end
