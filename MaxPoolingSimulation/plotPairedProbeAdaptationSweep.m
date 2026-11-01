function handles = plotPairedProbeAdaptationSweep(results)
% PLOTPAIREDPROBEADAPTATIONSWEEP  Gains, operating points, and winner anatomy.

T=results.conditionTable; W=results.winnerDirectionTable; cfg=results.config;
bases=unique(T.BaseDriftCoherence,'stable'); offsets=unique(T.OffsetDeg,'stable');
modes=["Independent","Yoked"]; colors=lines(numel(bases)); ocolors=lines(numel(offsets));
handles=gobjects(7,1);

handles(1)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation gains');
tl=tiledlayout(numel(offsets),2,'Padding','compact','TileSpacing','compact');
for oo=1:numel(offsets), for mm=1:2
    ax=nexttile((oo-1)*2+mm); hold(ax,'on');
    for bb=1:numel(bases)
        S=conditionRows(T,bases(bb),offsets(oo),modes(mm));
        plot(ax,S.AftereffectAt60,S.DriftGain,'-o','Color',colors(bb,:), ...
            'DisplayName',sprintf('Drift, base %.3g%%',bases(bb)));
        plot(ax,S.AftereffectAt60,S.ProbeGainPerStream,'--s','Color',colors(bb,:), ...
            'DisplayName',sprintf('Probe, base %.3g%%',bases(bb)));
    end
    yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
    xlabel(ax,'Aftereffect at 60% prestep (coherence-equivalent %)');
    ylabel(ax,'Recovered probit gain');
    title(ax,sprintf('%s, offset %.3g deg',modes(mm),offsets(oo)));
    if oo==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,'Recovered gains under drift-centered adaptation');

handles(2)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation operating points');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for mm=1:2, for rr=1:2
    ax=nexttile((rr-1)*2+mm); hold(ax,'on');
    for oo=1:numel(offsets), for bb=1:numel(bases)
        S=conditionRows(T,bases(bb),offsets(oo),modes(mm)); style='-';
        if bb==2, style='--'; end
        if rr==1, value=S.Threshold75; else, value=100*S.Performance; end
        plot(ax,S.AftereffectAt60,value,style,'Color',ocolors(oo,:), ...
            'DisplayName',sprintf('Offset %.3g, base %.3g%%',offsets(oo),bases(bb)));
    end, end
    if rr==1
        yline(ax,cfg.PairedConfig.FixedStepCoherence,':k','HandleVisibility','off');
        ylabel(ax,'75% threshold (coherence %)');
    else
        yline(ax,100*cfg.PairedConfig.TargetPerformance,':k','HandleVisibility','off');
        ylabel(ax,sprintf('Performance at %.3g%% (%%)', ...
            cfg.PairedConfig.FixedStepCoherence));
    end
    grid(ax,'on'); xlabel(ax,'Aftereffect at 60% prestep (coherence-equivalent %)');
    title(ax,modes(mm)); if rr==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,'Operating points with one fixed internal-noise SD');

handles(3)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation winner fractions');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases), for mm=1:2
    ax=nexttile((bb-1)*2+mm); hold(ax,'on');
    offset=nearest(unique(T.OffsetDeg),cfg.WinnerDirectionPlotOffsetDeg);
    S=conditionRows(T,bases(bb),offset,modes(mm));
    plot(ax,S.AftereffectAt60,S.ChangedProbeWinner,'-o','DisplayName','Changed: probe');
    plot(ax,S.AftereffectAt60,S.NoChangeProbeWinner,'--o','DisplayName','No-change: probe');
    plot(ax,S.AftereffectAt60,S.ChangedDriftWinner,'-s','DisplayName','Changed: drift');
    plot(ax,S.AftereffectAt60,S.NoChangeDriftWinner,'--s','DisplayName','No-change: drift');
    ylim(ax,[0 1]); grid(ax,'on'); ylabel(ax,'Winner fraction');
    xlabel(ax,'Aftereffect at 60% prestep (coherence-equivalent %)');
    title(ax,sprintf('%s, base %.3g%%, offset %.3g deg',modes(mm),bases(bb),offset));
    if bb==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,'Reallocation of drift- and probe-associated winners');

handles(4)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation fixed-winner totals');
tl=tiledlayout(numel(bases),2,'Padding','compact','TileSpacing','compact');
for bb=1:numel(bases), for mm=1:2
    ax=nexttile((bb-1)*2+mm); hold(ax,'on');
    offset=nearest(unique(W.OffsetDeg),cfg.WinnerDirectionPlotOffsetDeg);
    R=W(W.BaseDriftCoherence==bases(bb) & W.OffsetDeg==offset & W.Yoking==modes(mm),:);
    [~,first]=unique(R.AftereffectAt60,'stable'); S=sortrows(R(first,:),'AftereffectAt60');
    plot(ax,S.AftereffectAt60,S.DriftFixedWinnerTotal,'-o','DisplayName','Drift');
    plot(ax,S.AftereffectAt60,S.ProbeFixedWinnerTotalPerStream,'--s', ...
        'DisplayName','Probe per stream');
    yline(ax,0,':k','HandleVisibility','off'); grid(ax,'on');
    ylabel(ax,'Winner-weighted fixed-winner derivative');
    xlabel(ax,'Aftereffect at 60% prestep (coherence-equivalent %)');
    title(ax,sprintf('%s, base %.3g%%, offset %.3g deg',modes(mm),bases(bb),offset));
    if bb==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,'Signed contribution integrated over winner direction');

plotBase=nearest(unique(W.BaseDriftCoherence),cfg.WinnerDirectionPlotBase);
plotOffset=nearest(unique(W.OffsetDeg),cfg.WinnerDirectionPlotOffsetDeg);
available=unique(W.AftereffectAt60); targets=cfg.WinnerDirectionPlotAftereffects;
aes=arrayfun(@(x)nearest(available,x),targets); aes=unique(aes,'stable');
handles(5)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation winner directions');
tl=tiledlayout(numel(aes),2,'Padding','compact','TileSpacing','compact');
for aa=1:numel(aes), for mm=1:2
    ax=nexttile((aa-1)*2+mm); hold(ax,'on'); S=directionRows(W,plotBase,plotOffset,modes(mm),aes(aa));
    stairs(ax,S.WinnerDirectionDeg,S.ChangedWinnerFraction,'-','DisplayName','Changed patch');
    stairs(ax,S.WinnerDirectionDeg,S.NoChangeWinnerFraction,'--','DisplayName','No-change patch');
    stairs(ax,S.WinnerDirectionDeg,S.PatchAverageWinnerFraction,':','LineWidth',1.4, ...
        'DisplayName','Patch average');
    xlim(ax,[-180 180]); grid(ax,'on'); xlabel(ax,'Winner direction relative to drift (deg)');
    ylabel(ax,sprintf('Probability per %.3g-deg bin',cfg.WinnerDirectionBinWidthDeg));
    title(ax,sprintf('%s, aftereffect %.3g%%',modes(mm),aes(aa)));
    if aa==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,sprintf('Winner distributions: base %.3g%%, offset %.3g deg',plotBase,plotOffset));

handles(6)=figure('Visible',cfg.Visible,'Color','w','Name','Adaptation signed contributions');
tl=tiledlayout(numel(aes),2,'Padding','compact','TileSpacing','compact');
for aa=1:numel(aes), for mm=1:2
    ax=nexttile((aa-1)*2+mm); hold(ax,'on'); S=directionRows(W,plotBase,plotOffset,modes(mm),aes(aa));
    stairs(ax,S.WinnerDirectionDeg,S.DriftWinnerWeightedContribution,'-','DisplayName','Drift');
    stairs(ax,S.WinnerDirectionDeg,S.ProbeWinnerWeightedContributionPerStream,'--', ...
        'DisplayName','Probe per stream');
    yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[-180 180]); grid(ax,'on');
    xlabel(ax,'Winner direction relative to drift (deg)');
    ylabel(ax,'Winner probability x local sensitivity');
    title(ax,sprintf('%s, aftereffect %.3g%%: sum D=%.3f, P=%.3f',modes(mm), ...
        aes(aa),S.DriftFixedWinnerTotal(1),S.ProbeFixedWinnerTotalPerStream(1)));
    if aa==1 && mm==1, legend(ax,'Location','best'); end
end, end
title(tl,sprintf('Signed contributions: base %.3g%%, offset %.3g deg',plotBase,plotOffset));

% Direct empirical constraints: base-dependent threshold increment and the
% observed threshold-by-offset profiles.  The IDS 120-deg point is the primary
% 17-session anchor; the 170-deg point currently contains only 3 sessions.
B=results.baseContrastTable;
handles(7)=figure('Visible',cfg.Visible,'Color','w','Name','Empirical threshold constraints');
tl=tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
for mm=1:2
    ax=nexttile(mm); hold(ax,'on');
    for oo=1:numel(offsets)
        S=sortrows(B(B.OffsetDeg==offsets(oo) & B.Yoking==modes(mm),:), ...
            'AftereffectAt60');
        plot(ax,S.AftereffectAt60,S.ThresholdIncrease,'-o','Color',ocolors(oo,:), ...
            'DisplayName',sprintf('Offset %.3g deg',offsets(oo)));
    end
    yline(ax,cfg.ExpectedBaseThresholdIncrease,':k', ...
        sprintf('Provisional target %.3g%%',cfg.ExpectedBaseThresholdIncrease), ...
        'HandleVisibility','off');
    grid(ax,'on'); xlabel(ax,'Aftereffect at 60% prestep (coherence-equivalent %)');
    ylabel(ax,'Threshold increase: 40% base - 7.071% base'); title(ax,modes(mm));
    if mm==1, legend(ax,'Location','best'); end
end

selectedTargets=[0 4 7]; selectedAEs=arrayfun(@(x)nearest(available,x),selectedTargets);
selectedAEs=unique(selectedAEs,'stable'); aeColors=lines(numel(selectedAEs));
ax=nexttile(3); hold(ax,'on');
for aa=1:numel(selectedAEs)
    S=thresholdByOffset(T,bases(1),"Independent",selectedAEs(aa));
    plot(ax,S.OffsetDeg,S.Threshold75,'-o','Color',aeColors(aa,:), ...
        'DisplayName',sprintf('Model aftereffect %.3g%%',selectedAEs(aa)));
end
plot(ax,cfg.EmpiricalIDSOffsetsDeg,cfg.EmpiricalIDSThresholds,'kd', ...
    'MarkerFaceColor','k','DisplayName','IDS means');
for k=1:numel(cfg.EmpiricalIDSSessionCounts)
    if isfinite(cfg.EmpiricalIDSSessionCounts(k))
        text(ax,cfg.EmpiricalIDSOffsetsDeg(k),cfg.EmpiricalIDSThresholds(k), ...
            sprintf('  n=%g',cfg.EmpiricalIDSSessionCounts(k)), ...
            'VerticalAlignment','bottom');
    end
end
grid(ax,'on'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'75% threshold (coherence %)');
title(ax,'Low base, independent probes versus IDS'); legend(ax,'Location','best');

ax=nexttile(4); hold(ax,'on');
for aa=1:numel(selectedAEs)
    S=thresholdByOffset(T,bases(end),"Yoked",selectedAEs(aa));
    plot(ax,S.OffsetDeg,S.Threshold75,'-o','Color',aeColors(aa,:), ...
        'DisplayName',sprintf('Model aftereffect %.3g%%',selectedAEs(aa)));
end
yline(ax,cfg.EmpiricalIDRThresholdRange(1),':k','IDR mean range', ...
    'HandleVisibility','off');
yline(ax,cfg.EmpiricalIDRThresholdRange(2),':k','HandleVisibility','off');
grid(ax,'on'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'75% threshold (coherence %)');
title(ax,'High base, yoked probes versus IDR'); legend(ax,'Location','best');
title(tl,'Empirical threshold constraints on the adaptation model');
end

function S=conditionRows(T,base,offset,mode)
S=sortrows(T(T.BaseDriftCoherence==base & T.OffsetDeg==offset & T.Yoking==mode,:), ...
    'AftereffectAt60');
end
function S=directionRows(W,base,offset,mode,ae)
S=sortrows(W(W.BaseDriftCoherence==base & W.OffsetDeg==offset & W.Yoking==mode & ...
    W.AftereffectAt60==ae,:),'WinnerDirectionDeg');
end
function S=thresholdByOffset(T,base,mode,ae)
S=sortrows(T(T.BaseDriftCoherence==base & T.Yoking==mode & ...
    T.AftereffectAt60==ae,:),'OffsetDeg');
end
function value=nearest(values,target)
[~,k]=min(abs(values-target)); value=values(k);
end
