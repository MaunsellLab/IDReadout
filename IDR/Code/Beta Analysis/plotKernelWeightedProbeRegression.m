function fig = plotKernelWeightedProbeRegression(reg)
% plotKernelWeightedProbeRegression  Per-session coefficient diagnostics.

labels = {'DEC','INC','Combined'};
fields = {'dec','inc','combined'};
fig = figure('Color','w','Position',[100 100 1050 360]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for i = 1:3
  nexttile; hold on;
  F = reg.fitByStep.(fields{i});
  if isfield(F,'fitUsable') && F.fitUsable
    values = [F.betaPref F.betaProbe];
    errors = [F.betaPrefSE F.betaProbeSE];
    errorbar(1:2, values, errors, 'o', 'LineStyle','none', ...
      'MarkerFaceColor','auto','LineWidth',1.2);
    yline(0,':');
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'Preferred','Probe'});
    ylabel('\beta per % coherence');
    title(sprintf('%s, n=%d',labels{i},F.nTrials));
    txt = sprintf('probe/pref = %.3g\nSE = %.3g\ngrad = %.2g', ...
      F.betaRatio,F.betaRatioSE,F.gradientInfNorm);
    yl = ylim;
    text(0.58,yl(2)-0.08*range(yl),txt,'VerticalAlignment','top');
  else
    axis off;
    title(labels{i});
    if isfield(F,'message'), text(0.05,0.5,F.message,'Interpreter','none'); end
  end
  box off;
end

sessionID = '';
probeTag = '';
if isfield(reg.sessionProbeHeader,'sessionID'), sessionID = reg.sessionProbeHeader.sessionID; end
if isfield(reg.sessionProbeHeader,'probeTag'), probeTag = reg.sessionProbeHeader.probeTag; end
sgtitle(sprintf('%s %s kernel-weighted preferred/probe regression',sessionID,probeTag), ...
  'Interpreter','none');
end
