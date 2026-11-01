function fig = plotMTReadoutPNormGains(recovered,outputPath)
%PLOTMTREADOUTPNORMGAINS Compare observed and recovered IDS gain curves.
if isstruct(recovered), recovered=num2cell(recovered); end
assert(~isempty(recovered),'No recovered gains.');
x=cellfun(@(s)s.probeOffsetDeg,recovered); [x,order]=sort(x);
recovered=recovered(order); ps=recovered{1}.p;
for k=1:numel(recovered)
  assert(isequal(recovered{k}.p,ps),'Pooling grids differ.');
end
n=numel(x); obs=nan(n,2); ci=nan(n,2,2); ratio=nan(n,1); rci=nan(n,2);
for k=1:n
  o=recovered{k}.observed;
  obs(k,:)=o.driftProbe.gain(:)'; ci(k,:,:)=reshape(o.driftProbe.CI95,[1 2 2]);
  ratio(k)=o.normalized.pooledProbeOverDrift; rci(k,:)=o.normalized.pooledCI95;
end
fig=figure('Color','w'); tl=tiledlayout(2,1,'TileSpacing','compact');
colors=lines(numel(ps));
ax=nexttile(tl); hold(ax,'on');
errorbar(ax,x,ratio,ratio-rci(:,1),rci(:,2)-ratio,'ko', ...
  'MarkerFaceColor','k','DisplayName','Observed');
for j=1:numel(ps)
  y=cellfun(@(s)s.predicted{j}.normalized.pooledProbeOverDrift,recovered);
  plot(ax,x,y,'o-','Color',colors(j,:),'DisplayName',sprintf('Model p=%g',ps(j)));
end
yline(ax,0,'k:','HandleVisibility','off'); xlim(ax,[0 180]);
ylabel(ax,'Probe / drift gain'); xlabel(ax,'Probe offset (deg)');
grid(ax,'on'); box(ax,'off'); legend(ax,'Location','best');
ax=nexttile(tl); hold(ax,'on'); markers={'o','s'};
observedLabels=["Observed drift","Observed probe"];
for c=1:2
  errorbar(ax,x,obs(:,c),obs(:,c)-ci(:,c,1),ci(:,c,2)-obs(:,c), ...
    ['k' markers{c}],'DisplayName',observedLabels(c));
end
for j=1:numel(ps)
  y=nan(n,2);
  for k=1:n, y(k,:)=recovered{k}.predicted{j}.driftProbe.gain(:)'; end
  plot(ax,x,y(:,1),'--','Color',colors(j,:),'DisplayName',sprintf('Drift p=%g',ps(j)));
  plot(ax,x,y(:,2),'-','Color',colors(j,:),'DisplayName',sprintf('Probe p=%g',ps(j)));
end
yline(ax,0,'k:','HandleVisibility','off'); xlim(ax,[0 180]);
ylabel(ax,'Absolute noise gain'); xlabel(ax,'Probe offset (deg)');
grid(ax,'on'); box(ax,'off'); legend(ax,'Location','best');
title(tl,'MT population p-norm: observed and predicted IDS gains');
if nargin>=2 && strlength(string(outputPath))>0
  exportgraphics(fig,outputPath,'ContentType','vector');
end
end
