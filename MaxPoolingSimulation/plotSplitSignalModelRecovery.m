function figures = plotSplitSignalModelRecovery(results,cfg)
% PLOTSPLITSIGNALMODELRECOVERY  Threshold, recovery, and gain summaries.

nModels=numel(cfg.PoolingModels); colors=lines(nModels);
figures.thresholds=figure('Color','w','Visible',cfg.Visible); hold on;
for m=1:nModels
    plot(cfg.AllocationFractions,results.thresholds(m,:),'-o','LineWidth',1.5, ...
        'Color',colors(m,:),'DisplayName',modelLabel(cfg.PoolingModels{m}));
end
xlabel('Larger-stream allocation'); ylabel('75% threshold coherence (%)');
grid on; legend('Location','best'); title('Model-specific calibrated thresholds');

figures.recovery=figure('Color','w','Visible',cfg.Visible); hold on;
sessions=1:cfg.MaxRecoverySessions;
for m=1:nModels
    plot(sessions,results.recovery(m).diagonalProbability,'-o','LineWidth',1.5, ...
        'Color',colors(m,:),'DisplayName',modelLabel(cfg.PoolingModels{m}));
end
yline(cfg.TargetRecoveryProbability,':'); ylim([0 1]);
xlabel('Sessions'); ylabel('P(true model selected)'); grid on;
legend('Location','southeast'); title('Full refit model recovery');

figures.gains=figure('Color','w','Visible',cfg.Visible); tiledlayout(1,nModels);
for m=1:nModels
    nexttile; hold on; g=results.gains(m);
    errorbar(cfg.AllocationFractions,g.stream1Mean,1.96*g.stream1MeanSE, ...
        'o-','LineWidth',1.2,'DisplayName','Larger stream');
    errorbar(cfg.AllocationFractions,g.stream2Mean,1.96*g.stream2MeanSE, ...
        's-','LineWidth',1.2,'DisplayName','Smaller stream');
    yline(0,':'); xlabel('Larger-stream allocation'); ylabel('Mean frame gain');
    title(modelLabel(cfg.PoolingModels{m})); grid on;
    if m==1, legend('Location','best'); end
end
end

function label=modelLabel(name)
switch lower(name)
    case 'max', label='Max';
    case 'vector', label='Population vector magnitude';
    case 'energy', label='Population squared energy';
    otherwise, label=name;
end
end
