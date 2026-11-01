function groupedResult = runMaxPoolingDesignGroupFit( ...
    cachePath,resultInput,designProfilesInput,outputFolder,varargin)
%RUNMAXPOOLINGDESIGNGROUPFIT Fit the selected three-adaptation grouped model.

p=inputParser;
p.addParameter('DevelopmentDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('FinalDraws',256,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.parse(varargin{:});
if ~isfolder(outputFolder), mkdir(outputFolder); end
baseResult=loadInput(resultInput,'result');
designProfiles=loadInput(designProfilesInput,'designProfiles');
cache=maxpoolfit.loadCache(cachePath);
assert(isfield(designProfiles,'trialGroup') && ...
  numel(designProfiles.trialGroup)==height(cache.trials), ...
  'Saved design-profile labels do not match the current cache.');
trialGroup=string(designProfiles.trialGroup(:));
if isfield(designProfiles,'sourceDesign')
  sourceDesign=designProfiles.sourceDesign;
else
  [sourceDesign,auditGroup]=maxpoolfit.driftRegimeMetadata(cache);
  assert(isequal(auditGroup,trialGroup), ...
    'Current session classification differs from the saved design profiles.');
end
assert(isfield(designProfiles,'config'), ...
  'Design profiles do not contain their fitted model configuration.');
cfg=designProfiles.config;
if isfield(cfg,'coherenceResponseMode') && ...
    string(cfg.coherenceResponseMode)=="random-dot-redistribution"
  cfg.modelVersion="DenseHardMax-v20-redistribution-design-group-adaptation";
else
  cfg.modelVersion="DenseHardMax-v18-design-group-adaptation";
end

P=designProfiles.groupProfile;
required=["IDQ-fixed" "IDQ-randomized" "IDS" "IDR"];
assert(all(ismember(required,unique(string(P.analysisGroup)))), ...
  'Design profiles do not contain all four required groups.');
aFixed=bestGroupAdaptation(P,"IDQ-fixed");
aIDR=bestGroupAdaptation(P,"IDR");
aShared=bestSharedAdaptation(P,"IDQ-randomized","IDS");
analysisGroup=["IDQ-fixed";"IDQ-randomized";"IDS";"IDR"];
adaptation=[aFixed;aShared;aShared;aIDR];
adaptationByGroup=table(analysisGroup,adaptation);

developmentFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'TrialMask',baseResult.finalFit.trialMask, ...
  'NDraws',p.Results.DevelopmentDraws);
finalFit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
  adaptationByGroup,'TrialMask',baseResult.finalFit.trialMask, ...
  'NDraws',p.Results.FinalDraws);
prediction=maxpoolfit.predictDesignGroupModel(cache,finalFit,cfg, ...
  'NDraws',p.Results.FinalDraws);
groupedResult=struct('developmentFit',developmentFit,'finalFit',finalFit, ...
  'prediction',prediction,'sourceDesign',sourceDesign, ...
  'designProfileSource',designProfiles,'cachePath',string(cachePath), ...
  'createdAt',datetime('now'));
groupedResult.adaptationSelection=struct( ...
  'IDQFixed',aFixed,'IDQRandomizedAndIDS',aShared,'IDR',aIDR, ...
  'source',sprintf('%d-draw design-group profile grid', ...
    designProfiles.nNoiseDraws));

result=groupedResult; %#ok<NASGU>
save(fullfile(outputFolder,'MaxPoolingDesignGroupFitResults.mat'), ...
  'result','-v7.3');
writetable(finalFit.groupFit,fullfile(outputFolder, ...
  'MaxPoolingDesignGroupFit_Summary.csv'));
fig=figure('Color','w','Position',[100 100 900 430]);
labels=categorical(finalFit.groupFit.analysisGroup, ...
  finalFit.groupFit.analysisGroup,'Ordinal',true);
yyaxis left; bar(labels,finalFit.groupFit.noiseSDPC);
ylabel('Fitted internal-noise SD (% coherence)');
yyaxis right; plot(labels,finalFit.groupFit.adaptation,'ko-', ...
  'MarkerFaceColor','k','LineWidth',1.2);
ylabel('Fixed adaptation, a'); ylim([0 1]); grid on; box off;
title('Grouped hard-max fit: selected adaptation and refitted noise');
exportgraphics(fig,fullfile(outputFolder, ...
  'MaxPoolingDesignGroupFit_Parameters.pdf'),'ContentType','vector');
close(fig);
end

function value=loadInput(input,name)
if isstruct(input), value=input; return; end
S=load(input,name);
assert(isfield(S,name),'MAT file does not contain variable %s.',name);
value=S.(name);
end

function a=bestGroupAdaptation(P,group)
Q=P(string(P.analysisGroup)==group,:);
[~,row]=max(Q.logLikelihood); a=Q.adaptation(row);
end

function a=bestSharedAdaptation(P,group1,group2)
A=P(string(P.analysisGroup)==group1,:);
B=P(string(P.analysisGroup)==group2,:);
[grid,ia,ib]=intersect(A.adaptation,B.adaptation,'stable');
assert(~isempty(grid),'The two group profiles have no common adaptation grid.');
[~,row]=max(A.logLikelihood(ia)+B.logLikelihood(ib));
a=grid(row);
end
