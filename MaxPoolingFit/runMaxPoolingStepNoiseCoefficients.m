function stepCoefficients = runMaxPoolingStepNoiseCoefficients( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGSTEPNOISECOEFFICIENTS Diagnose RDK gains at exact step sizes.

p=inputParser;
p.addParameter('MinimumTrials',100, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('ReliableSessionCount',5, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
cache=maxpoolfit.loadCache(cachePath);
if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit') && isfield(result,'prediction'), ...
  'A completed primary fit and prediction are required.');

stepCoefficients=maxpoolfit.noiseCoefficientsByStep( ...
  cache,result.prediction,'TrialMask',result.finalFit.trialMask, ...
  'MinimumTrials',p.Results.MinimumTrials, ...
  'ReliableSessionCount',p.Results.ReliableSessionCount);
stepCoefficients.modelVersion=result.finalFit.modelVersion;
stepCoefficients.cachePath=string(cachePath);
stepCoefficients.winFractionAgreement= ...
  summarizeWinFractions(stepCoefficients.summary);
save(fullfile(outputFolder,'MaxPoolingStepNoiseCoefficients.mat'), ...
  'stepCoefficients','-v7.3');
writetable(stepCoefficients.summary,fullfile(outputFolder, ...
  'MaxPoolingStepNoiseCoefficients.csv'));
writetable(stepCoefficients.winFractionAgreement,fullfile(outputFolder, ...
  'MaxPoolingGainImpliedWinFractionAgreement.csv'));
plotSource(stepCoefficients,"observed",outputFolder);
plotSource(stepCoefficients,"primary domain noise",outputFolder);
plotGainWinFractions(stepCoefficients,outputFolder);
end

function A=summarizeWinFractions(S)
sources=["observed" "primary domain noise"];
components=["Drift" "Probe"];
A=table();
for s=1:numel(sources)
  for k=1:numel(components)
    prefix=lower(components(k));
    flagName=prefix+"GainSignsConsistent";
    valueName=prefix+"GainImpliedWinFraction";
    q=string(S.source)==sources(s) & S.inferenceReliable & S.(flagName);
    G=S(q,:); weights=double(G.nTrials);
    source=sources(s); component=components(k); nCells=height(G);
    nTrials=sum(weights);
    if isempty(G)
      weightedAccuracy=nan; weightedGainImpliedWinFraction=nan;
      weightedBias=nan; weightedRMSE=nan;
    else
      weightedAccuracy=sum(weights.*G.accuracy)./sum(weights);
      weightedGainImpliedWinFraction= ...
        sum(weights.*G.(valueName))./sum(weights);
      differences=G.(valueName)-G.accuracy;
      weightedBias=sum(weights.*differences)./sum(weights);
      weightedRMSE=sqrt(sum(weights.*differences.^2)./sum(weights));
    end
    row=table(source,component,nCells,nTrials,weightedAccuracy, ...
      weightedGainImpliedWinFraction,weightedBias,weightedRMSE);
    if isempty(A), A=row; else, A=[A;row]; end %#ok<AGROW>
  end
end
end

function plotGainWinFractions(C,outputFolder)
S=C.summary;
groups=unique(S(:,{'domain','animal'}),'rows','stable');
colors=lines(height(groups));
sources=["observed" "primary domain noise"];
sourceTitles=["Observed choices" "Primary model"];
components=["Drift" "Probe"];
fig=figure('Color','w','Position',[100 100 1050 820]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for s=1:numel(sources)
  for k=1:numel(components)
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    plot(ax,[0.5 0.9],[0.5 0.9],'k--','LineWidth',1, ...
      'HandleVisibility','off');
    prefix=lower(components(k));
    flagName=prefix+"GainSignsConsistent";
    valueName=prefix+"GainImpliedWinFraction";
    seName=valueName+"SE";
    nShown=0;
    for g=1:height(groups)
      q=string(S.source)==sources(s) & S.inferenceReliable & ...
        S.(flagName) & string(S.domain)==string(groups.domain(g)) & ...
        string(S.animal)==string(groups.animal(g));
      G=S(q,:); if isempty(G), continue; end
      error=1.96.*G.(seName); error(~isfinite(error))=0;
      errorbar(ax,G.accuracy,G.(valueName),error,'o', ...
        'LineStyle','none','Color',colors(g,:), ...
        'MarkerFaceColor',colors(g,:),'MarkerSize',5,'CapSize',3, ...
        'DisplayName',sprintf('%s-%s',groups.domain(g),groups.animal(g)));
      values=G.(valueName);
      for j=1:height(G)
        text(ax,G.accuracy(j)+0.004,values(j), ...
          sprintf('%g/%g',G.probeOffsetDeg(j),G.stepMagnitudePC(j)), ...
          'FontSize',6,'Color',colors(g,:),'Clipping','on');
      end
      nShown=nShown+height(G);
    end
    xlim(ax,[0.5 0.9]); ylim(ax,[0 1]);
    xlabel(ax,'Accuracy');
    ylabel(ax,'Gain-implied change-RDK win fraction');
    title(ax,sprintf('%s - %s noise (n=%d)', ...
      sourceTitles(s),lower(components(k)),nShown));
    if s==1 && k==1, legend(ax,'Location','southwest'); end
  end
end
title(tl,['Gain-implied win fraction versus accuracy' newline ...
  'point labels: probe offset / exact step magnitude (% coherence)']);
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingGainImpliedWinFraction.pdf'),'ContentType','vector');
close(fig);
end

function plotSource(C,source,outputFolder)
S=C.summary(string(C.summary.source)==source,:);
groups=unique(S(:,{'domain','animal'}),'rows','stable');
fig=figure('Color','w','Position',[100 100 1100 max(700,230*height(groups))]);
tl=tiledlayout(height(groups),2,'TileSpacing','compact','Padding','compact');
for g=1:height(groups)
  q=string(S.domain)==string(groups.domain(g)) & ...
    string(S.animal)==string(groups.animal(g));
  G=S(q,:); offsets=unique(G.probeOffsetDeg,'sorted'); colors=lines(numel(offsets));
  ax=nexttile(tl); plotNoise(ax,G,offsets,colors,"Drift");
  ylabel(ax,sprintf('%s - %s',groups.domain(g),groups.animal(g)));
  if g==1
    title(ax,'Drift noise: solid changed; dashed -unchanged');
    legend(ax,'Location','best','NumColumns',min(4,numel(offsets)));
  end
  ax=nexttile(tl); plotNoise(ax,G,offsets,colors,"Probe");
  if g==1, title(ax,'Probe noise: solid changed; dashed -unchanged'); end
end
xlabel(tl,'Step magnitude (% coherence)');
ylabel(tl,'Choice coefficient per % coherence');
title(tl,sprintf('%s coefficients within exact step conditions',upperFirst(source)));
if source=="observed"
  fileName='MaxPoolingObservedNoiseCoefficientsByStep.pdf';
else
  fileName='MaxPoolingPredictedNoiseCoefficientsByStep.pdf';
end
exportgraphics(fig,fullfile(outputFolder,fileName),'ContentType','vector');
close(fig);
end

function plotNoise(ax,G,offsets,colors,kind)
hold(ax,'on'); grid(ax,'on'); box(ax,'off'); yline(ax,0,'k:','HandleVisibility','off');
changedName="changed"+kind+"CoefficientPerPC";
unchangedName="unchanged"+kind+"CoefficientPerPC";
changedSE="changed"+kind+"CoefficientSE";
unchangedSE="unchanged"+kind+"CoefficientSE";
for k=1:numel(offsets)
  R=sortrows(G(G.probeOffsetDeg==offsets(k),:),'stepMagnitudePC');
  changedError=1.96.*R.(changedSE); changedError(~isfinite(changedError))=0;
  unchangedError=1.96.*R.(unchangedSE); unchangedError(~isfinite(unchangedError))=0;
  errorbar(ax,R.stepMagnitudePC,R.(changedName),changedError, ...
    'o-','Color',colors(k,:),'MarkerFaceColor',colors(k,:), ...
    'LineWidth',1.15,'CapSize',3,'DisplayName',sprintf('%g deg',offsets(k)));
  errorbar(ax,R.stepMagnitudePC,-R.(unchangedName),unchangedError, ...
    's--','Color',colors(k,:),'MarkerFaceColor','w', ...
    'LineWidth',1.05,'CapSize',3,'HandleVisibility','off');
end
end

function value=upperFirst(value)
value=string(value); value=upper(extractBefore(value,2))+extractAfter(value,1);
end
