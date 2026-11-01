function gainProfile2D = runMaxPoolingAdaptationGain2DProfile( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGADAPTATIONGAIN2DPROFILE Sparse joint profile of a and b.

p=inputParser;
p.addParameter('AdaptationGrid',[0 0.3 0.56 0.8 1], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&& ...
  all(x>=0)&&all(x<=1));
p.addParameter('GainGrid',[0 0.1 0.25 0.5 1], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('ReferenceAdaptation',[], ...
  @(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1));
p.addParameter('ReferencePC',40, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
if isstruct(resultInput)
  result=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  result=S.result;
end
assert(isfield(result,'finalFit')&&isfinite(result.finalFit.adaptation), ...
  'A completed primary fit is required.');
cache=maxpoolfit.loadCache(cachePath);
cfg=result.finalFit.config;
if isempty(p.Results.ReferenceAdaptation)
  referenceAdaptation=result.finalFit.adaptation;
else
  referenceAdaptation=double(p.Results.ReferenceAdaptation);
end
adaptationGrid=unique([double(p.Results.AdaptationGrid(:)); ...
  referenceAdaptation],'sorted');
gainGrid=unique([double(p.Results.GainGrid(:));0],'sorted');
cfg.adaptation=referenceAdaptation;
cfg=maxpoolfit.configureAdaptationExponent(cfg,0,p.Results.ReferencePC);
cfg=maxpoolfit.configureAdaptationGain(cfg,0,p.Results.ReferencePC);
gainProfile2D=maxpoolfit.profileAdaptationGain(cache,cfg, ...
  'TrialMask',result.finalFit.trialMask,'NDraws',p.Results.NDraws, ...
  'Adaptation',referenceAdaptation,'AdaptationGrid',adaptationGrid, ...
  'GainGrid',gainGrid,'ReferencePC',p.Results.ReferencePC);
stem='MaxPoolingIncrementAdaptationGain2DProfile';
save(fullfile(outputFolder,[stem '.mat']),'gainProfile2D','-v7.3');

P=gainProfile2D.profile;
exportTable=P(:,{'adaptation','adaptationGainStrength','logLikelihood', ...
  'deltaLogLFromBest','deltaLogLFromLinear'});
for id=1:numel(gainProfile2D.domains)
  suffix=matlab.lang.makeValidName(gainProfile2D.domains(id));
  exportTable.("noiseSDPC_"+suffix)=P.noiseSDPC(:,id);
  exportTable.("logLikelihood_"+suffix)=P.domainLogLikelihood(:,id);
  exportTable.("evalCount_"+suffix)=P.evalCount(:,id);
end
writetable(exportTable,fullfile(outputFolder,[stem '.csv']));

Z=nan(numel(adaptationGrid),numel(gainGrid));
for ia=1:numel(adaptationGrid)
  for ib=1:numel(gainGrid)
    q=abs(P.adaptation-adaptationGrid(ia))<1e-12 & ...
      abs(P.adaptationGainStrength-gainGrid(ib))<1e-12;
    Z(ia,ib)=P.deltaLogLFromBest(q);
  end
end
fig=figure('Color','w','Position',[100 100 760 560]);
imagesc(gainGrid,adaptationGrid,Z); axis xy; colorbar;
hold on; contour(gainGrid,adaptationGrid,Z,[-20 -10 -5 -2 0], ...
  'k','ShowText','on');
plot(P.adaptation(gainProfile2D.bestRow), ...
  P.adaptationGainStrength(gainProfile2D.bestRow),'wp', ...
  'MarkerFaceColor','w','MarkerSize',12);
xlabel('Increment gain strength, b'); ylabel('Subtractive adaptation, a');
title(sprintf('Joint adaptation profile; P_0 = %.1f%%',p.Results.ReferencePC));
exportgraphics(fig,fullfile(outputFolder,[stem '.pdf']),'ContentType','vector');
close(fig);
end
