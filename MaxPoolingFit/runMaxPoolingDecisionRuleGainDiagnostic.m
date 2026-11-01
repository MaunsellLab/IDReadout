function decisionRuleGain = runMaxPoolingDecisionRuleGainDiagnostic( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGDECISIONRULEGAINDIAGNOSTIC Test explicit prestep subtraction.

if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
result=loadResult(resultInput);
decisionRuleGain=maxpoolfit.decisionRuleGainDiagnostic( ...
  cache,result,varargin{:});
decisionRuleGain.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingDecisionRuleGainDiagnostic.mat'), ...
  'decisionRuleGain','-v7.3');
writetable(decisionRuleGain.gainSummary,fullfile(outputFolder, ...
  'MaxPoolingDecisionRuleGainDiagnostic.csv'));
writetable(decisionRuleGain.modelComparison,fullfile(outputFolder, ...
  'MaxPoolingDecisionRuleModelComparison.csv'));
writetable(decisionRuleGain.domainComparison,fullfile(outputFolder, ...
  'MaxPoolingDecisionRuleDomainComparison.csv'));
plotDiagnostic(decisionRuleGain,outputFolder);
end

function result=loadResult(input)
if isstruct(input), result=input; return; end
S=load(input,'result');
assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
result=S.result;
end

function plotDiagnostic(D,outputFolder)
S=D.gainSummary; labels=strrep(S.source,'-',' ');
changedDrift=S.changedDriftCoefficientPerPC;
unchangedDrift=-S.unchangedDriftCoefficientPerPC;
changedProbe=S.changedProbeCoefficientPerPC;
unchangedProbe=-S.unchangedProbeCoefficientPerPC;
fig=figure('Color','w','Position',[100 100 1100 470]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); b=bar(ax,[changedDrift unchangedDrift]); grid(ax,'on'); box(ax,'off');
ax.XTick=1:height(S); ax.XTickLabel=labels; ax.XTickLabelRotation=18;
ylabel(ax,'Choice coefficient per % coherence'); title(ax,'Drift');
legend(b,{'Changed RDK','-Unchanged RDK'},'Location','best');
ax=nexttile(tl); bar(ax,[changedProbe unchangedProbe]); grid(ax,'on'); box(ax,'off');
ax.XTick=1:height(S); ax.XTickLabel=labels; ax.XTickLabelRotation=18;
ylabel(ax,'Choice coefficient per % coherence'); title(ax,'Probe');
Q=D.selection;
title(tl,sprintf('Decision-rule gain diagnostic: %s - %s; %g deg; %g%% step', ...
  Q.domain,Q.animal,Q.probeOffsetDeg,Q.stepMagnitudePC));
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingDecisionRuleGainDiagnostic.pdf'),'ContentType','vector');
close(fig);
end
