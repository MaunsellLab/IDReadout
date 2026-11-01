function opponentAxis=runMaxPoolingOpponentAxisDecomposition( ...
    cachePath,jointFitInput,outputFolder,varargin)
%RUNMAXPOOLINGOPPONENTAXISDECOMPOSITION Diagnose paired opponent probes.

if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
jointFit=loadJointFit(jointFitInput);
opponentAxis=maxpoolfit.opponentAxisDecomposition(cache,jointFit,varargin{:});
opponentAxis.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingOpponentAxisDecomposition.mat'), ...
  'opponentAxis','-v7.3');
writetable(opponentAxis.summary,fullfile(outputFolder, ...
  'MaxPoolingOpponentAxisDecomposition_Summary.csv'));
writetable(opponentAxis.trialSummary,fullfile(outputFolder, ...
  'MaxPoolingOpponentAxisDecomposition_Trials.csv'));
writetable(opponentAxis.gainSummary,fullfile(outputFolder, ...
  'MaxPoolingOpponentAxisDecomposition_Gains.csv'));
plotSummary(opponentAxis,outputFolder);
end

function value=loadJointFit(input)
if isstruct(input), value=input; return; end
S=load(input);
if isfield(S,'jointGainFit'), value=S.jointGainFit;
else
  names=fieldnames(S); assert(numel(names)==1, ...
    'Joint-fit MAT file must contain one result variable.');
  value=S.(names{1});
end
assert(isfield(value,'gainOptimalFit')&&isfield(value,'choiceOptimalFit'), ...
  'Input is not a completed joint gain fit.');
end

function plotSummary(R,outputFolder)
S=R.summary; fits=unique(S.fitSelection,'stable');
fig=figure('Color','w','Position',[70 70 1150 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
metrics={"meanDriftMinusOpponentResponse","internalProbeWinnerFraction", ...
  "probeNoiseChangesMechanismFraction","predictedChangedProbeGain"};
titles={"Deterministic drift-opponent gap", ...
  "Probe controls within-RDK maximum", ...
  "Experimental probe noise changes winner", ...
  "Probe influence after between-RDK comparison"};
ylabs={"Response difference","Fraction of internal-noise draws", ...
  "Fraction of internal-noise draws","Predicted changed-side probe gain"};
for m=1:numel(metrics)
  ax=nexttile(tl); hold(ax,'on');
  configs=unique(S.configuration,'stable'); x=1:numel(configs);
  for f=1:numel(fits)
    Q=S(S.fitSelection==fits(f),:);
    [~,loc]=ismember(configs,Q.configuration);
    style='o-'; if f==2, style='s--'; end
    plot(ax,x,Q.(metrics{m})(loc),style,'LineWidth',1.4, ...
      'DisplayName',fits(f));
  end
  xticks(ax,x); xticklabels(ax,strrep(configs,' ','\newline')); grid(ax,'on');
  yline(ax,0,'k:','HandleVisibility','off'); ylabel(ax,ylabs{m});
  title(ax,titles{m}); if m==1, legend(ax,'Location','best'); end
end
title(tl,sprintf('Opponent-axis decomposition: IDS %g deg (%d trials)', ...
  R.probeOffsetDeg,R.nTrials));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingOpponentAxisDecomposition.pdf'),'ContentType','vector'); close(fig);
end
