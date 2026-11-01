function widthFit=runMaxPoolingDirectionalWidthGainFit( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGDIRECTIONALWIDTHGAINFIT Fit effective width to changed gains.

p=inputParser;
p.addParameter('DirectionalWidthGridDeg',[10 15 20 25 30 37.5 45 60],@isnumeric);
p.addParameter('TargetGroups',["IDQ-randomized" "IDS"], ...
  @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MinimumGainSessions',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath); result=loadInput(resultInput,'result');
assert(isfield(result,'finalFit')&&isfield(result.finalFit,'adaptationByGroup'), ...
  'A completed design-group fit is required.');
widthFit=maxpoolfit.profileDirectionalWidthByGain(cache,result.finalFit, ...
  'DirectionalWidthGridDeg',p.Results.DirectionalWidthGridDeg, ...
  'TargetGroups',p.Results.TargetGroups,'NDraws',p.Results.NDraws, ...
  'MinimumGainSessions',p.Results.MinimumGainSessions);
widthFit.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingDirectionalWidthGainFit.mat'), ...
  'widthFit','-v7.3');
writetable(widthFit.overallProfile,fullfile(outputFolder, ...
  'MaxPoolingDirectionalWidthGainFit_Overall.csv'));
writetable(widthFit.groupProfile,fullfile(outputFolder, ...
  'MaxPoolingDirectionalWidthGainFit_Groups.csv'));
writetable(widthFit.gainProfile,fullfile(outputFolder, ...
  'MaxPoolingDirectionalWidthGainFit_Gains.csv'));
plotObjectives(widthFit,outputFolder);
plotGainCurves(widthFit,outputFolder);
end

function value=loadInput(input,name)
if isstruct(input), value=input; return; end
S=load(input,name); assert(isfield(S,name), ...
  'MAT file does not contain variable %s.',name); value=S.(name);
end

function plotObjectives(P,outputFolder)
O=P.overallProfile; fig=figure('Color','w','Position',[80 80 1050 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); plot(ax,O.directionalWidthDeg,O.changedGainRMSE,'ko-', ...
  'MarkerFaceColor','k','LineWidth',1.3); grid(ax,'on');
xline(ax,P.gainOptimalWidthDeg,'b--','Gain optimum');
xlabel(ax,'Effective directional sigma (deg)'); ylabel(ax,'Changed-gain RMSE');
title(ax,'Primary gain objective');
ax=nexttile(tl); plot(ax,O.directionalWidthDeg, ...
  O.deltaTargetLogLFromChoiceBest,'ko-','MarkerFaceColor','k','LineWidth',1.3);
grid(ax,'on'); xline(ax,P.choiceOptimalWidthDeg,'r--','Choice optimum');
xlabel(ax,'Effective directional sigma (deg)'); ylabel(ax,'Delta log likelihood');
title(ax,'Target-group choice likelihood');
ax=nexttile(tl); hold(ax,'on');
plot(ax,O.directionalWidthDeg,O.changedDriftGainRMSE,'o-', ...
  'DisplayName','drift','LineWidth',1.2);
plot(ax,O.directionalWidthDeg,O.changedProbeGainRMSE,'o-', ...
  'DisplayName','probe','LineWidth',1.2);
grid(ax,'on'); legend(ax,'Location','best'); xlabel(ax,'Effective directional sigma (deg)');
ylabel(ax,'Gain RMSE'); title(ax,'Components of gain error');
ax=nexttile(tl); plot(ax,O.directionalWidthDeg,O.meanTargetNoiseSDPC,'ko-', ...
  'MarkerFaceColor','k','LineWidth',1.3); grid(ax,'on');
xlabel(ax,'Effective directional sigma (deg)'); ylabel(ax,'Mean fitted noise SD (% coherence)');
title(ax,'Choice-performance nuisance parameter');
title(tl,sprintf('Gain-targeted directional-width profile; %d draws',P.nNoiseDraws));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingDirectionalWidthGainFit_Objectives.pdf'),'ContentType','vector');
close(fig);
end

function plotGainCurves(P,outputFolder)
G=P.gainProfile; widths=[P.gainOptimalWidthDeg P.choiceOptimalWidthDeg];
labels=["gain-optimal" "choice-optimal"];
if widths(1)==widths(2), widths=widths(1); labels="joint optimum"; end
groups=P.targetGroups; fig=figure('Color','w','Position',[60 60 1200 780]);
tl=tiledlayout(2,numel(groups),'TileSpacing','compact','Padding','compact');
for component=["changedDriftGain" "changedProbeGain"]
  for g=1:numel(groups)
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    O=G(string(G.domain)==groups(g)&string(G.source)=="observed"& ...
      G.directionalWidthDeg==P.directionalWidthGridDeg(1),:);
    animals=unique(string(O.animal),'stable'); colors=lines(numel(animals));
    for a=1:numel(animals)
      q=string(O.animal)==animals(a); [x,order]=sort(O.probeOffsetDeg(q));
      y=O.(component)(q); y=y(order); se=O.(component+"SE")(q); se=se(order);
      se(~isfinite(se))=0;
      errorbar(ax,x,y,se,'o','Color',colors(a,:), ...
        'MarkerFaceColor',colors(a,:),'DisplayName',animals(a)+" observed");
      for w=1:numel(widths)
        M=G(string(G.domain)==groups(g)& ...
          string(G.source)=="directional-width model"& ...
          G.directionalWidthDeg==widths(w)&string(G.animal)==animals(a),:);
        [mx,mo]=sort(M.probeOffsetDeg); my=M.(component); my=my(mo);
        style='-'; if w==2, style='--'; end
        plot(ax,mx,my,style,'Color',colors(a,:),'LineWidth',1.5, ...
          'DisplayName',animals(a)+" "+labels(w));
      end
    end
    yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)');
    ylabel(ax,strrep(component,'changed','changed '));
    title(ax,groups(g)+" - "+strrep(component,"changed","")+" ");
    if component=="changedDriftGain"&&g==1, legend(ax,'Location','best'); end
  end
end
title(tl,sprintf('Observed and fitted changed-side gains; gain optimum = %g deg', ...
  P.gainOptimalWidthDeg));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingDirectionalWidthGainFit_Curves.pdf'),'ContentType','vector');
close(fig);
end
