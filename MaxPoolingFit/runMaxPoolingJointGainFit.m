function jointGainFit=runMaxPoolingJointGainFit( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGJOINTGAINFIT Joint gain fit of width, adaptation, correlation.

if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
baseFit=resolveDesignGroupFit(resultInput,outputFolder);
checkpointPath=fullfile(outputFolder,'MaxPoolingJointGainFit_Checkpoint.mat');
jointGainFit=maxpoolfit.optimizeGainMechanism(cache,baseFit, ...
  'CheckpointPath',checkpointPath,varargin{:});
jointGainFit.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingJointGainFit.mat'), ...
  'jointGainFit','-v7.3');
writetable(jointGainFit.evaluation,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Evaluations.csv'));
writetable(jointGainFit.groupEvaluation,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Groups.csv'));
writetable(jointGainFit.gainEvaluation,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Gains.csv'));
writetable(jointGainFit.paretoEvaluation,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Pareto.csv'));
plotObjectives(jointGainFit,outputFolder);
plotGainCurves(jointGainFit,outputFolder);
if isfile(checkpointPath), delete(checkpointPath); end
end

function fit=resolveDesignGroupFit(input,outputFolder)
value=[];
if isstruct(input)
  value=input;
elseif isfile(input)
  S=load(input);
  names=fieldnames(S);
  assert(numel(names)==1, ...
    'Input MAT file must contain one result variable.');
  value=S.(names{1});
end
fit=extractFit(value);
if ~isempty(fit), return; end

% The ordinary calling pattern often leaves the original primary result in
% the workspace while the required grouped result is saved in outputFolder.
groupedPath=fullfile(outputFolder,'MaxPoolingDesignGroupFitResults.mat');
if isfile(groupedPath)
  S=load(groupedPath,'result');
  if isfield(S,'result'), fit=extractFit(S.result); end
end
assert(~isempty(fit),[ ...
  'A completed design-group fit is required. Pass the output of ' ...
  'runMaxPoolingDesignGroupFit, or its saved ' ...
  'MaxPoolingDesignGroupFitResults.mat file.']);
end

function fit=extractFit(value)
fit=[];
if ~isstruct(value), return; end
if isfield(value,'adaptationByGroup')&&isfield(value,'trialGroup')
  fit=value;
elseif isfield(value,'finalFit')&& ...
    isfield(value.finalFit,'adaptationByGroup')
  fit=value.finalFit;
elseif isfield(value,'gainOptimalFit')&& ...
    isfield(value.gainOptimalFit,'adaptationByGroup')
  fit=value.gainOptimalFit;
elseif isfield(value,'choiceOptimalFit')&& ...
    isfield(value.choiceOptimalFit,'adaptationByGroup')
  fit=value.choiceOptimalFit;
end
end

function plotObjectives(P,outputFolder)
E=P.evaluation; fig=figure('Color','w','Position',[80 80 1120 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); scatter(ax,E.changedGainRMSE,E.targetLogLikelihood,35, ...
  E.directionalWidthDeg,'filled'); hold(ax,'on');
Q=E.isPareto; scatter(ax,E.changedGainRMSE(Q),E.targetLogLikelihood(Q), ...
  65,'ko','LineWidth',1.3); grid(ax,'on'); colorbar(ax);
xlabel(ax,'Changed-gain RMSE'); ylabel(ax,'Target log likelihood');
title(ax,'Gain-choice tradeoff (color: width)');
ax=nexttile(tl); scatter3(ax,E.directionalWidthDeg,E.adaptation, ...
  E.correlationLengthDeg,45,E.changedGainRMSE,'filled'); grid(ax,'on');
xlabel(ax,'Width (deg)'); ylabel(ax,'Adaptation'); zlabel(ax,'Correlation (deg)');
colorbar(ax); title(ax,'Search geometry (color: gain RMSE)');
ax=nexttile(tl); scatter(ax,E.directionalWidthDeg,E.changedGainRMSE,35, ...
  E.adaptation,'filled'); grid(ax,'on'); colorbar(ax);
xlabel(ax,'Width (deg)'); ylabel(ax,'Changed-gain RMSE');
title(ax,'Width and adaptation');
ax=nexttile(tl); scatter(ax,E.correlationLengthDeg,E.changedGainRMSE,35, ...
  E.meanTargetNoiseSDPC,'filled'); grid(ax,'on'); colorbar(ax);
xlabel(ax,'Correlation length (deg)'); ylabel(ax,'Changed-gain RMSE');
title(ax,'Correlation and fitted noise SD');
title(tl,sprintf('Joint gain fit; %d search draws, %d final draws', ...
  P.nNoiseDraws,P.finalNNoiseDraws));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Objectives.pdf'),'ContentType','vector'); close(fig);
end

function plotGainCurves(P,outputFolder)
G=[finalGains(P.gainOptimalComparison,"gain-optimal"); ...
  finalGains(P.choiceOptimalComparison,"choice-optimal")];
groups=P.targetGroups; fig=figure('Color','w','Position',[60 60 1200 780]);
tl=tiledlayout(2,numel(groups),'TileSpacing','compact','Padding','compact');
for component=["changedDriftGain" "changedProbeGain"]
  for g=1:numel(groups)
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    O=G(string(G.domain)==groups(g)&string(G.source)=="observed",:);
    O=unique(O,'rows','stable'); animals=unique(string(O.animal),'stable');
    colors=lines(numel(animals));
    for a=1:numel(animals)
      q=string(O.animal)==animals(a); [x,order]=sort(O.probeOffsetDeg(q));
      y=O.(component)(q); y=y(order); se=O.(component+"SE")(q); se=se(order);
      errorbar(ax,x,y,se,'o','Color',colors(a,:), ...
        'MarkerFaceColor',colors(a,:),'DisplayName',animals(a)+" observed");
      for source=["gain-optimal" "choice-optimal"]
        M=G(string(G.domain)==groups(g)&string(G.source)==source& ...
          string(G.animal)==animals(a),:);
        [mx,mo]=sort(M.probeOffsetDeg); my=M.(component); my=my(mo);
        style='-'; if source=="choice-optimal", style='--'; end
        plot(ax,mx,my,style,'Color',colors(a,:),'LineWidth',1.5, ...
          'DisplayName',animals(a)+" "+source);
      end
    end
    yline(ax,0,'k:','HandleVisibility','off'); xlabel(ax,'Probe offset (deg)');
    ylabel(ax,strrep(component,'changed','changed '));
    title(ax,groups(g)+" - "+strrep(component,"changed","")+" ");
    if component=="changedDriftGain"&&g==1, legend(ax,'Location','best'); end
  end
end
title(tl,'Observed gains and joint-fit predictions');
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingJointGainFit_Curves.pdf'),'ContentType','vector'); close(fig);
end

function G=finalGains(comparison,label)
G=comparison.summary;
G=G(ismember(string(G.source),["observed" "primary domain noise"]),:);
G.source(string(G.source)=="primary domain noise")=label;
end
