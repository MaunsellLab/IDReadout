function diagnostics = runMaxPoolingMechanismDiagnostics( ...
    cachePath,resultInput,performanceInput,outputFolder,varargin)
%RUNMAXPOOLINGMECHANISMDIAGNOSTICS Diagnose near-opponent max competition.

p=inputParser;
p.addParameter('ConditionRows',[],@(x)isnumeric(x)&&isvector(x));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MaxTrialsPerCondition',1000, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.addParameter('NoiseGridPC',[0 1 2 3 4 5 6 8 10 12 15 20 30 40],@isnumeric);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
result=loadStructInput(resultInput,'result');
performanceMatched=loadStructInput(performanceInput,'performanceMatched');
assert(string(result.finalFit.modelVersion)== ...
  string(performanceMatched.fit.baseModelVersion), ...
  'Primary and performance-matched inputs use different model versions.');
assert(abs(result.finalFit.adaptation-performanceMatched.fit.adaptation)<1e-12, ...
  'Primary and performance-matched fits use different adaptation values.');
[primaryQ,primaryB,primaryP0]=adaptationWeighting(result.finalFit);
[matchedQ,matchedB,matchedP0]=adaptationWeighting(performanceMatched.fit);
assert(abs(primaryQ-matchedQ)<1e-12 && abs(primaryB-matchedB)<1e-12 && ...
  abs(primaryP0-matchedP0)<1e-12, ...
  'Primary and performance-matched fits use different adaptation weighting.');
G=performanceMatched.fit.groupFit;
if isempty(p.Results.ConditionRows)
  rows=defaultConditionRows(G);
else
  rows=p.Results.ConditionRows(:);
end
diagnostics=maxpoolfit.mechanismDiagnostics( ...
  cache,result.finalFit,performanceMatched.fit,rows, ...
  'NDraws',p.Results.NDraws, ...
  'MaxTrialsPerCondition',p.Results.MaxTrialsPerCondition, ...
  'NoiseGridPC',p.Results.NoiseGridPC);
save(fullfile(outputFolder,'MaxPoolingMechanismDiagnostics.mat'), ...
  'diagnostics','-v7.3');
writetable(diagnostics.conditionSummary,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_Conditions.csv'));
writetable(diagnostics.noiseSweep,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_NoiseSweep.csv'));
writetable(diagnostics.winnerHistogram,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_WinnerDirections.csv'));
plotNoiseSweep(diagnostics,outputFolder);
plotMechanismFractions(diagnostics,outputFolder);
plotWinnerDirections(diagnostics,outputFolder);
end

function value=loadStructInput(input,name)
if isstruct(input)
  value=input;
else
  S=load(input,name);
  assert(isfield(S,name),'MAT file does not contain variable %s.',name);
  value=S.(name);
end
end

function [q,b,P0]=adaptationWeighting(fit)
q=0; b=0; P0=40;
if isfield(fit,'adaptationExponent'), q=fit.adaptationExponent;
elseif isfield(fit,'config')&&isfield(fit.config,'adaptationExponent')
  q=fit.config.adaptationExponent;
end
if isfield(fit,'adaptationGainStrength'), b=fit.adaptationGainStrength;
elseif isfield(fit,'config')&&isfield(fit.config,'adaptationGainStrength')
  b=fit.config.adaptationGainStrength;
end
if isfield(fit,'adaptationReferencePC'), P0=fit.adaptationReferencePC;
elseif isfield(fit,'config')&&isfield(fit.config,'adaptationReferencePC')
  P0=fit.config.adaptationReferencePC;
end
end

function rows=defaultConditionRows(G)
target=find(G.nTrials>=500 & ...
  ((string(G.domain)=="IDS" & G.probeOffsetDeg>=160) | ...
   (string(G.domain)=="IDR" & G.probeOffsetDeg>=170)));
assert(~isempty(target),'No default near-opponent conditions were found.');
rows=target(:);
for k=1:numel(target)
  q=string(G.domain)==string(G.domain(target(k))) & ...
    string(G.animal)==string(G.animal(target(k))) & ...
    abs(G.probeOffsetDeg-45)<1e-5 & G.nTrials>=500;
  candidates=find(q);
  if ~isempty(candidates)
    [~,best]=min(abs(G.stepMagnitudePC(candidates)-G.stepMagnitudePC(target(k))));
    rows(end+1,1)=candidates(best); %#ok<AGROW>
  end
end
rows=unique(rows,'stable');
end

function labels=conditionLabels(S)
labels=string(S.domain)+" - "+string(S.animal)+" - "+ ...
  string(round(S.probeOffsetDeg))+" deg - "+ ...
  string(round(S.stepMagnitudePC,2))+"%";
end

function plotNoiseSweep(D,outputFolder)
S=D.conditionSummary; W=D.noiseSweep; labels=conditionLabels(S);
n=height(S); nCols=min(3,n); nRows=ceil(n/nCols);
fig=figure('Color','w','Position',[100 100 1100 300*nRows]);
tl=tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');
for k=1:n
  ax=nexttile(tl); q=conditionRowMask(W,S(k,:));
  plot(ax,W.noiseSDPC(q),100.*W.predictedAccuracy(q),'o-','LineWidth',1.2); hold(ax,'on');
  yline(ax,100*S.observedAccuracy(k),'k--','Observed');
  xline(ax,S.primaryNoiseSDPC(k),':','Primary');
  xline(ax,S.performanceNoiseSDPC(k),'--','Matched');
  grid(ax,'on'); box(ax,'off');
  xlabel(ax,'Noise SD (% coherence)'); ylabel(ax,'Predicted accuracy (%)');
  title(ax,labels(k));
end
title(tl,'Attainable accuracy as internal noise varies');
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_NoiseSweep.pdf'),'ContentType','vector');
close(fig);
end

function plotMechanismFractions(D,outputFolder)
S=D.conditionSummary; labels=conditionLabels(S);
x=categorical(labels,labels,'Ordinal',true);
fig=figure('Color','w','Position',[100 100 1100 500]);
bar(x,100.*[S.changedDriftWinFraction S.unchangedDriftWinFraction ...
  S.bothRDKsProbeWinFraction]);
ylabel('Trial-draws (%)'); grid on; box off;
legend({'Drift region wins: changed RDK','Drift region wins: unchanged RDK', ...
  'Probe region wins: both RDKs'},'Location','best');
title('Within-RDK competition at performance-matched noise');
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_Competition.pdf'),'ContentType','vector');
close(fig);
end

function plotWinnerDirections(D,outputFolder)
S=D.conditionSummary; H=D.winnerHistogram; labels=conditionLabels(S);
n=height(S); nCols=min(3,n); nRows=ceil(n/nCols);
fig=figure('Color','w','Position',[100 100 1100 300*nRows]);
tl=tiledlayout(nRows,nCols,'TileSpacing','compact','Padding','compact');
for k=1:n
  ax=nexttile(tl); q=conditionRowMask(H,S(k,:));
  plot(ax,H.winnerBinCenterDeg(q),100.*H.changedWinnerFraction(q), ...
    '-','LineWidth',1.3); hold(ax,'on');
  plot(ax,H.winnerBinCenterDeg(q),100.*H.unchangedWinnerFraction(q), ...
    '-','LineWidth',1.3);
  xline(ax,0,'k:'); xlim(ax,[-180 180]); grid(ax,'on'); box(ax,'off');
  xlabel(ax,'Winning direction relative to drift (deg)');
  ylabel(ax,'Winners (%)'); title(ax,labels(k));
end
lgd=legend(ax,{'Changed RDK','Unchanged RDK'}, ...
  'Orientation','horizontal');
lgd.Layout.Tile='south';
title(tl,'Within-RDK winning mechanism directions');
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingMechanismDiagnostics_WinnerDirections.pdf'),'ContentType','vector');
close(fig);
end

function q=conditionRowMask(T,row)
q=string(T.domain)==string(row.domain) & ...
  string(T.animal)==string(row.animal) & ...
  abs(T.probeOffsetDeg-row.probeOffsetDeg)<1e-5 & ...
  abs(T.stepMagnitudePC-row.stepMagnitudePC)<1e-5;
end
