function gainComparison = runMaxPoolingGainComparison( ...
    cachePath,resultInput,performanceInput,outputFolder,varargin)
%RUNMAXPOOLINGGAINCOMPARISON Compare observed and model choice-noise gains.

p=inputParser;
p.addParameter('TrialMask',[],@(x)isempty(x)||islogical(x)||isnumeric(x));
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
result=loadStructInput(resultInput,'result');
performanceMatched=loadStructInput(performanceInput,'performanceMatched');
assert(string(result.finalFit.modelVersion)== ...
  string(performanceMatched.fit.baseModelVersion), ...
  'Primary and performance-matched results use different model versions.');
isGrouped=isfield(result.finalFit,'adaptationByGroup');
if isGrouped
  assert(isfield(performanceMatched.fit,'adaptationByGroup'), ...
    'Performance-matched result lacks the design-group adaptations.');
  assert(equalAdaptationTables(result.finalFit.adaptationByGroup, ...
    performanceMatched.fit.adaptationByGroup), ...
    'Primary and performance-matched results use different group adaptations.');
else
  assert(abs(result.finalFit.adaptation-performanceMatched.fit.adaptation)<1e-12, ...
    'Primary and performance-matched results use different adaptation values.');
end
[primaryQ,primaryB,primaryP0]=adaptationWeighting(result.finalFit);
[matchedQ,matchedB,matchedP0]=adaptationWeighting(performanceMatched.fit);
assert(abs(primaryQ-matchedQ)<1e-12 && abs(primaryB-matchedB)<1e-12 && ...
  abs(primaryP0-matchedP0)<1e-12, ...
  'Primary and performance-matched results use different adaptation weighting.');
if isempty(p.Results.TrialMask), mask=result.finalFit.trialMask;
else, mask=p.Results.TrialMask; end

workCache=cache;
if isGrouped
  assert(numel(result.finalFit.trialGroup)==height(cache.trials), ...
    'Grouped fit labels do not match the cache.');
  workCache.trials.domain=result.finalFit.trialGroup;
end
gainComparison=maxpoolfit.gainComparison(workCache,result.prediction, ...
  performanceMatched.prediction,'TrialMask',mask);
primarySource="primary domain noise";
if isGrouped
  primarySource="primary design-group noise";
  gainComparison.summary.source( ...
    string(gainComparison.summary.source)=="primary domain noise")=primarySource;
  gainComparison.sources( ...
    string(gainComparison.sources)=="primary domain noise")=primarySource;
end
gainComparison.modelVersion=result.finalFit.modelVersion;
gainComparison.adaptation=result.finalFit.adaptation;
if isGrouped
  gainComparison.adaptationByGroup=result.finalFit.adaptationByGroup;
end
gainComparison.adaptationExponent=primaryQ;
gainComparison.adaptationGainStrength=primaryB;
gainComparison.adaptationReferencePC=primaryP0;
gainComparison.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingGainComparison.mat'), ...
  'gainComparison','-v7.3');
writetable(gainComparison.summary,fullfile(outputFolder, ...
  'MaxPoolingGainComparison.csv'));
writetable(gainComparison.referencePsychometric,fullfile(outputFolder, ...
  'MaxPoolingGainComparison_ReferencePsychometric.csv'));
plotGainSide(gainComparison,"changed",outputFolder);
plotGainSide(gainComparison,"unchanged",outputFolder);
plotChangedAcrossDomains(gainComparison,primarySource, ...
  outputFolder,'MaxPoolingGainComparison_ChangedAcrossDomains_Primary.pdf');
plotChangedAcrossDomains(gainComparison,"performance matched", ...
  outputFolder,'MaxPoolingGainComparison_ChangedAcrossDomains_PerformanceMatched.pdf');
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

function tf=equalAdaptationTables(A,B)
required={'analysisGroup','adaptation'};
if ~all(ismember(required,A.Properties.VariableNames)) || ...
    ~all(ismember(required,B.Properties.VariableNames))
  tf=false; return;
end
A=sortrows(A(:,required),'analysisGroup');
B=sortrows(B(:,required),'analysisGroup');
tf=height(A)==height(B) && ...
  isequal(string(A.analysisGroup),string(B.analysisGroup)) && ...
  all(abs(A.adaptation-B.adaptation)<1e-12);
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

function plotGainSide(C,side,outputFolder)
S=C.summary; groups=unique(S(:,{'domain','animal'}),'rows','stable');
nGroups=height(groups);
fig=figure('Color','w','Position',[100 100 1050 max(650,230*nGroups)]);
tl=tiledlayout(nGroups,2,'TileSpacing','compact','Padding','compact');
for g=1:nGroups
  q=string(S.domain)==string(groups.domain(g)) & ...
    string(S.animal)==string(groups.animal(g));
  G=S(q,:);
  ax=nexttile(tl); plotGain(ax,G,side+"DriftGain",side+"DriftGainSE");
  ylabel(ax,sprintf('%s - %s',groups.domain(g),groups.animal(g)));
  if g==1
    title(ax,'Drift noise');
    legend(ax,cellstr(C.sources),'Location','best');
  end
  ax=nexttile(tl); plotGain(ax,G,side+"ProbeGain",side+"ProbeGainSE");
  if g==1, title(ax,'Probe noise'); end
end
xlabel(tl,'Probe offset (degrees)');
ylabel(tl,'Equivalent step-coherence gain');
title(tl,sprintf('%s-RDK observed versus predicted gains', ...
  upperFirst(side)));
fileName="MaxPoolingGainComparison_"+upperFirst(side)+"RDK.pdf";
exportgraphics(fig,fullfile(outputFolder,fileName),'ContentType','vector');
close(fig);
end

function plotGain(ax,G,gainName,seName)
hold(ax,'on'); grid(ax,'on'); box(ax,'off'); yline(ax,0,'k:','HandleVisibility','off');
sourceOrder=unique(string(G.source),'stable').';
colors=[0.1 0.1 0.1; 0.1 0.4 0.8; 0.9 0.4 0.1];
markers={'o','s','d'};
for s=1:3
  q=string(G.source)==sourceOrder(s); R=sortrows(G(q,:),'probeOffsetDeg');
  x=R.probeOffsetDeg; y=R.(gainName);
  if s==1
    plot(ax,x,y,[markers{s} '-'], ...
      'Color',colors(s,:),'MarkerFaceColor',colors(s,:), ...
      'LineWidth',1.1,'DisplayName',sourceOrder(s));
    reliable=R.nSessions>=3; observedSE=R.(seName);
    errorbar(ax,x(reliable),y(reliable),1.96.*observedSE(reliable), ...
      'LineStyle','none','Color',colors(s,:),'CapSize',4, ...
      'HandleVisibility','off');
    markPreliminary(ax,R,x,y,colors(s,:));
  else
    plot(ax,x,y,[markers{s} '-'],'Color',colors(s,:), ...
      'MarkerFaceColor',colors(s,:),'LineWidth',1.2, ...
      'DisplayName',sourceOrder(s));
  end
end
offsets=unique(G.probeOffsetDeg); pad=max(5,0.04.*max(offsets));
if numel(offsets)==1, xlim(ax,[offsets-10 offsets+10]);
else, xlim(ax,[min(offsets)-pad max(offsets)+pad]); end
end

function value=upperFirst(value)
value=string(value); value=upper(extractBefore(value,2))+extractAfter(value,1);
end

function markPreliminary(ax,R,x,y,color)
q=R.nSessions<3;
if ~any(q), return; end
plot(ax,x(q),y(q),'x','Color',color,'LineWidth',1.6, ...
  'MarkerSize',8,'HandleVisibility','off');
for k=find(q(:).')
  text(ax,x(k),y(k),sprintf('  %d session%s',R.nSessions(k), ...
    pluralS(R.nSessions(k))),'Color',color,'FontSize',8, ...
    'VerticalAlignment','bottom','Clipping','on');
end
end

function s=pluralS(n)
if n==1, s=''; else, s='s'; end
end

function plotChangedAcrossDomains(C,modelSource,outputFolder,fileName)
% Compact progress plot for the primary changed-side scientific target.
S=C.summary; animals=unique(string(S.animal),'stable');
domains=unique(string(S.domain),'stable');
colors=lines(numel(domains));
markerCycle={'o','s','d','^','v','>','<','p'};
variables=["changedDriftGain" "changedProbeGain" "changedProbeOverDrift"];
titles=["Drift gain" "Probe gain" "Probe / drift gain"];
fig=figure('Color','w','Position',[70 70 1280 max(520,390*numel(animals))]);
tl=tiledlayout(numel(animals),3,'TileSpacing','compact','Padding','compact');
for ia=1:numel(animals)
  for iv=1:numel(variables)
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    yline(ax,0,'k:','HandleVisibility','off');
    for id=1:numel(domains)
      marker=markerCycle{1+mod(id-1,numel(markerCycle))};
      base=string(S.animal)==animals(ia) & string(S.domain)==domains(id);
      O=sortrows(S(base & string(S.source)=="observed",:),'probeOffsetDeg');
      M=sortrows(S(base & string(S.source)==modelSource,:),'probeOffsetDeg');
      if isempty(O)||isempty(M), continue; end
      assert(isequal(O.probeOffsetDeg,M.probeOffsetDeg), ...
        'Observed and fitted offset grids differ for %s - %s.', ...
        domains(id),animals(ia));
      if iv<=2
        seName=replace(variables(iv),"Gain","GainSE");
        plot(ax,O.probeOffsetDeg,O.(variables(iv)), ...
          'LineStyle','-','Marker',marker,'Color',colors(id,:), ...
          'MarkerFaceColor',colors(id,:),'LineWidth',1.15, ...
          'DisplayName',domains(id)+" observed");
        reliable=O.nSessions>=3; observedValue=O.(variables(iv));
        observedSE=O.(seName);
        errorbar(ax,O.probeOffsetDeg(reliable),observedValue(reliable), ...
          1.96.*observedSE(reliable),'LineStyle','none', ...
          'Color',colors(id,:),'CapSize',4,'HandleVisibility','off');
      else
        plot(ax,O.probeOffsetDeg,O.(variables(iv)),'LineStyle','-', ...
          'Marker',marker,'Color',colors(id,:), ...
          'MarkerFaceColor',colors(id,:),'LineWidth',1.15, ...
          'DisplayName',domains(id)+" observed");
      end
      markPreliminary(ax,O,O.probeOffsetDeg,O.(variables(iv)),colors(id,:));
      plot(ax,M.probeOffsetDeg,M.(variables(iv)),'LineStyle','--', ...
        'Marker',marker,'Color',colors(id,:), ...
        'MarkerFaceColor','w','LineWidth',1.35, ...
        'DisplayName',domains(id)+" fit");
    end
    xlim(ax,[-5 185]); xticks(ax,[0 25 45 90 120 135 170 180]);
    if ia==1, title(ax,titles(iv)); end
    if ia==numel(animals), xlabel(ax,'Probe offset (degrees)'); end
    if iv==1, ylabel(ax,animals(ia)+"    Equivalent gain"); end
    if ia==1 && iv==3, legend(ax,'Location','best','NumColumns',2); end
  end
end
modelLabel=strrep(char(modelSource),'-',' ');
title(tl,sprintf('Changed-RDK gains across domains: observed versus %s',modelLabel));
exportgraphics(fig,fullfile(outputFolder,fileName),'ContentType','vector');
close(fig);
end
