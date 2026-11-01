function temporalGain = runMaxPoolingTemporalGainDiagnostics( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGTEMPORALGAINDIAGNOSTICS Five-bin RDK gain diagnostic.
%
% The default selection is the well-sampled IDR-Meetz 25-degree condition.

p=inputParser;
p.addParameter('Domain',"IDR",@(x)isstring(x)||ischar(x));
p.addParameter('Animal',"Meetz",@(x)isstring(x)||ischar(x));
p.addParameter('ProbeOffsetDeg',25,@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
p.addParameter('StepMagnitudePC',[],@(x)isempty(x)||(isnumeric(x)&&isvector(x)));
p.addParameter('BinFrames',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&mod(x,1)==0);
p.addParameter('NBootstrap',200,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&mod(x,1)==0);
p.addParameter('MinimumTrials',500,@(x)isnumeric(x)&&isscalar(x)&&x>=20&&mod(x,1)==0);
p.addParameter('MinimumSessions',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('RandomSeed',271828,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
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
assert(isfield(result,'finalFit')&&isfield(result,'prediction'), ...
  'A completed primary fit and prediction are required.');

temporalGain=maxpoolfit.temporalGainDiagnostics( ...
  cache,result.prediction,'TrialMask',result.finalFit.trialMask, ...
  'Domain',p.Results.Domain,'Animal',p.Results.Animal, ...
  'ProbeOffsetDeg',p.Results.ProbeOffsetDeg, ...
  'StepMagnitudePC',p.Results.StepMagnitudePC, ...
  'BinFrames',p.Results.BinFrames,'NBootstrap',p.Results.NBootstrap, ...
  'MinimumTrials',p.Results.MinimumTrials,'MinimumSessions',p.Results.MinimumSessions, ...
  'RandomSeed',p.Results.RandomSeed);
temporalGain.modelVersion=result.finalFit.modelVersion;
temporalGain.cachePath=string(cachePath);
save(fullfile(outputFolder,'MaxPoolingTemporalGainDiagnostics.mat'), ...
  'temporalGain','-v7.3');
writetable(temporalGain.profile,fullfile(outputFolder,'MaxPoolingTemporalGainProfiles.csv'));
writetable(temporalGain.shapeTests,fullfile(outputFolder,'MaxPoolingTemporalGainShapeTests.csv'));
plotProfiles(temporalGain,false,outputFolder);
plotProfiles(temporalGain,true,outputFolder);
end

function plotProfiles(D,normalized,outputFolder)
sources=["observed" "primary domain noise"]; components=["drift" "probe"];
fig=figure('Color','w','Position',[100 100 1050 760]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for s=1:2
  for c=1:2
    ax=nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'off');
    G=D.profile(string(D.profile.source)==sources(s) & ...
      string(D.profile.component)==components(c),:);
    C=sortrows(G(string(G.rdkStatus)=="changed",:),'timeBin');
    U=sortrows(G(string(G.rdkStatus)=="unchanged",:),'timeBin');
    x=(C.timeStartMS+C.timeEndMS)./2;
    if normalized
      changed=C.normalizedCoefficient; unchanged=U.normalizedCoefficient;
      cLow=max(0,changed-C.bootstrapLower95); cHigh=max(0,C.bootstrapUpper95-changed);
      uLow=max(0,unchanged-U.bootstrapLower95); uHigh=max(0,U.bootstrapUpper95-unchanged);
      errorbar(ax,x,changed,cLow,cHigh,'o-','LineWidth',1.4, ...
        'MarkerFaceColor',[0.1 0.35 0.75],'Color',[0.1 0.35 0.75], ...
        'DisplayName','Changed RDK');
      errorbar(ax,x,unchanged,uLow,uHigh,'s--','LineWidth',1.4, ...
        'MarkerFaceColor','w','Color',[0.75 0.2 0.25], ...
        'DisplayName','-Unchanged RDK');
      ylabel(ax,'Fraction of whole-step coefficient');
    else
      errorbar(ax,x,C.signedCoefficientPerPC,1.96.*C.coefficientSE, ...
        'o-','LineWidth',1.4,'MarkerFaceColor',[0.1 0.35 0.75], ...
        'Color',[0.1 0.35 0.75],'DisplayName','Changed RDK');
      errorbar(ax,x,U.signedCoefficientPerPC,1.96.*U.coefficientSE, ...
        's--','LineWidth',1.4,'MarkerFaceColor','w', ...
        'Color',[0.75 0.2 0.25],'DisplayName','-Unchanged RDK');
      yline(ax,0,'k:','HandleVisibility','off');
      ylabel(ax,'Choice coefficient per % coherence');
    end
    xlabel(ax,'Time in 250-ms step (ms)'); xlim(ax,[0 250]);
    R=D.shapeTests(string(D.shapeTests.source)==sources(s) & ...
      string(D.shapeTests.component)==components(c),:);
    title(ax,sprintf('%s - %s; shape p = %.3g', ...
      sourceTitle(sources(s)),components(c),R.shapePValue));
    if s==1&&c==1, legend(ax,'Location','best'); end
  end
end
selectionTitle=sprintf('%s - %s; probe offset %s deg; n=%d, sessions=%d', ...
  strjoin(D.selection.domain,','),strjoin(D.selection.animal,','), ...
  strjoin(string(D.selection.probeOffsetDeg),','),D.nTrials,D.nSessions);
if normalized
  title(tl,"Normalized temporal gain profiles"+newline+selectionTitle);
  fileName='MaxPoolingTemporalGainNormalizedProfiles.pdf';
else
  title(tl,"Five-bin temporal gain coefficients"+newline+selectionTitle);
  fileName='MaxPoolingTemporalGainProfiles.pdf';
end
exportgraphics(fig,fullfile(outputFolder,fileName),'ContentType','vector'); close(fig);
end

function titleText=sourceTitle(source)
if source=="observed", titleText="Observed"; else, titleText="Primary model"; end
end
