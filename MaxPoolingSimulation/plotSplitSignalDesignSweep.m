function figures = plotSplitSignalDesignSweep(summary,cfg)
% PLOTSPLITSIGNALDESIGNSWEEP  Primary performance, gain, and power plots.

rows=summary.primary; separations=unique([rows.separationDeg]);
allocations=unique([rows.allocation]); steps=unique([rows.totalStep]);
colors=lines(numel(allocations));

figures.performance=figure('Color','w','Visible',cfg.Visible); tiledlayout(1,numel(separations));
for s=1:numel(separations)
    nexttile; hold on;
    for a=1:numel(allocations)
        values=nan(size(steps));
        for j=1:numel(steps)
            hit=find([rows.separationDeg]==separations(s) & ...
                [rows.allocation]==allocations(a) & [rows.totalStep]==steps(j),1);
            if ~isempty(hit), values(j)=rows(hit).performance; end
        end
        plot(steps,values,'-o','Color',colors(a,:),'LineWidth',1.5, ...
            'DisplayName',sprintf('%.1f:%.1f',allocations(a),1-allocations(a)));
    end
    yline(cfg.TargetPerformance,':'); xlabel('Total step coherence (%)');
    ylabel('Proportion correct'); title(sprintf('%g deg separation',separations(s)));
    grid on; legend('Location','southeast');
end

figures.gains=figure('Color','w','Visible',cfg.Visible); tiledlayout(1,numel(separations));
for s=1:numel(separations)
    nexttile; hold on;
    for a=1:numel(allocations)
        candidates=rows([rows.separationDeg]==separations(s) & ...
            [rows.allocation]==allocations(a));
        [~,i]=min(abs([candidates.performance]-cfg.TargetPerformance)); r=candidates(i);
        plot(allocations(a),r.winnerGain,'o','Color',[0.1 0.45 0.8], ...
            'MarkerFaceColor',[0.1 0.45 0.8]);
        plot(allocations(a),r.loserGain,'s','Color',[0.8 0.25 0.2], ...
            'MarkerFaceColor',[0.8 0.25 0.2]);
    end
    yline(0,':'); xlabel('Larger-stream allocation'); ylabel('Probit noise gain');
    title(sprintf('%g deg: nearest 75%% correct',separations(s))); grid on;
end

figures.power=figure('Color','w','Visible',cfg.Visible); hold on;
powerField=sprintf('sep%d',round(cfg.PrimaryPowerSeparationDeg));
names={'maxMinusLinear','maxMinusVector','maxMinusEnergy'};
labels={'Linear','Vector magnitude','Squared energy'};
for k=1:numel(names)
    p=summary.power.(powerField).(names{k}); plot(p.sessions,p.probability,'-o','LineWidth',1.5, ...
        'DisplayName',labels{k});
end
yline(cfg.TargetPower,':'); xlabel('Sessions'); ylabel('P(DeltaLL criterion reached)');
ylim([0 1]); grid on; legend('Location','southeast');
title(sprintf('Primary design: %g deg separation',cfg.PrimaryPowerSeparationDeg));
end
