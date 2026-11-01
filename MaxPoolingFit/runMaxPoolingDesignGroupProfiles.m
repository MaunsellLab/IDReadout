function designProfiles = runMaxPoolingDesignGroupProfiles( ...
    cachePath,resultInput,outputFolder,varargin)
%RUNMAXPOOLINGDESIGNGROUPPROFILES Profile four design-defined trial groups.
%
% Groups: IDR, IDS, IDQ-randomized, and IDQ-fixed. IDQ classification uses
% sessionHeader.randomBaseDir with cached drift directions as an audit.

p=inputParser;
p.addParameter('AdaptationGrid',0:0.05:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('ReadHeaders',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('CoherenceResponseMode',"random-dot-redistribution", ...
  @(x)isstring(x)||ischar(x));
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
if isstruct(resultInput)
  fitResult=resultInput;
else
  S=load(resultInput,'result');
  assert(isfield(S,'result'),'Result MAT file does not contain variable result.');
  fitResult=S.result;
end
assert(isfield(fitResult,'finalFit')&&isfield(fitResult.finalFit,'trialMask'), ...
  'A completed primary fit result is required.');

cache=maxpoolfit.loadCache(cachePath);
cfg=fitResult.finalFit.config;
cfg=maxpoolfit.configureCoherenceResponse(cfg,p.Results.CoherenceResponseMode);
[sourceDesign,trialGroup]=maxpoolfit.driftRegimeMetadata(cache, ...
  'ReadHeaders',p.Results.ReadHeaders,'HeaderField',"randomBaseDir");
stem='MaxPoolingDesignGroupProfiles';
writetable(sourceDesign,fullfile(outputFolder,[stem '_SessionDesign.csv']));
disp(groupsummary(sourceDesign,'analysisGroup','sum','nCachedTrials'));
designProfiles=maxpoolfit.profileDesignGroups(cache,cfg,trialGroup, ...
  'TrialMask',fitResult.finalFit.trialMask, ...
  'NDraws',p.Results.NDraws, ...
  'AdaptationGrid',p.Results.AdaptationGrid);
designProfiles.sourceDesign=sourceDesign;

save(fullfile(outputFolder,[stem '.mat']),'designProfiles','-v7.3');
writetable(designProfiles.groupSummary, ...
  fullfile(outputFolder,[stem '_Summary.csv']));
writetable(designProfiles.groupProfile, ...
  fullfile(outputFolder,[stem '_Likelihood.csv']));
writetable(designProfiles.modelComparison, ...
  fullfile(outputFolder,[stem '_Comparisons.csv']));

fig=figure('Color','w','Position',[100 100 1050 450]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
colors=lines(numel(designProfiles.groupOrder));
for g=1:numel(designProfiles.groupOrder)
  Q=designProfiles.groupProfile.analysisGroup==designProfiles.groupOrder(g);
  plot(designProfiles.groupProfile.adaptation(Q), ...
    designProfiles.groupProfile.deltaLogLFromBest(Q),'o-', ...
    'Color',colors(g,:),'LineWidth',1.4, ...
    'DisplayName',designProfiles.groupOrder(g));
end
yline(0,'k:','HandleVisibility','off'); grid on; box off;
xlabel('Adaptation, a'); ylabel('\Delta log likelihood from group maximum');
title('Independent design-group profiles'); legend('Location','best');

nexttile;
S=designProfiles.groupSummary;
scatter(S.bestAdaptation,S.noiseSDPC,70,colors,'filled'); hold on;
for g=1:height(S)
  text(S.bestAdaptation(g),S.noiseSDPC(g),"  "+S.analysisGroup(g), ...
    'Color',colors(g,:),'VerticalAlignment','middle');
end
grid on; box off; xlabel('Best adaptation, a');
ylabel('Fitted internal-noise SD (% coherence)');
title('Parameter-space comparison');
exportgraphics(fig,fullfile(outputFolder,[stem '.pdf']),'ContentType','vector');
close(fig);
end
