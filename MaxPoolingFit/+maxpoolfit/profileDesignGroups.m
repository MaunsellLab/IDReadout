function result = profileDesignGroups(cache,cfg,trialGroup,varargin)
%PROFILEDESIGNGROUPS Profile IDR, IDS, randomized IDQ, and fixed IDQ.
%
% The four independent profiles are followed by nested descriptive
% comparisons within IDQ and between randomized IDQ and IDS.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('AdaptationGrid',0:0.05:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
trialGroup=string(trialGroup(:));
assert(numel(trialGroup)==height(cache.trials), ...
  'trialGroup must contain one label per cached trial.');
if isempty(p.Results.TrialMask)
  primaryMask=cfg.primarySelector(cache.trials);
else
  primaryMask=false(height(cache.trials),1);
  if islogical(p.Results.TrialMask)
    assert(numel(p.Results.TrialMask)==height(cache.trials), ...
      'Logical TrialMask must have one element per cached trial.');
    primaryMask=p.Results.TrialMask(:);
  else
    primaryMask(p.Results.TrialMask)=true;
  end
end

% Fit IDQ subsets first because their comparison determines whether the
% randomization manipulation matters before IDS is brought into the test.
expectedGroupOrder=["IDQ-fixed";"IDQ-randomized";"IDS";"IDR"];
present=false(size(expectedGroupOrder));
for g=1:numel(expectedGroupOrder)
  present(g)=any(primaryMask & trialGroup==expectedGroupOrder(g));
end
groupOrder=expectedGroupOrder(present);
assert(~isempty(groupOrder),'No recognized design-group trials were selected.');
if any(~present)
  warning('maxpoolfit:profileDesignGroups:MissingGroups', ...
    'Skipping design group(s) with no primary trials: %s.', ...
    strjoin(expectedGroupOrder(~present),', '));
end
profiles=cell(numel(groupOrder),1);
summaryRows=cell(0,9);
profileRows=cell(0,7);
for g=1:numel(groupOrder)
  groupMask=primaryMask & trialGroup==groupOrder(g);
  profiles{g}=maxpoolfit.profileAdaptation(cache,cfg, ...
    'TrialMask',groupMask,'NDraws',p.Results.NDraws, ...
    'AdaptationGrid',p.Results.AdaptationGrid);
  P=profiles{g}.profile;
  assert(size(P.noiseSDPC,2)==1, ...
    'Each design group must contain exactly one original domain.');
  domain=unique(string(cache.trials.domain(groupMask)),'stable');
  animal=string(cache.trials.animal(groupMask));
  session=string(cache.trials.sessionID(groupMask));
  sourceCount=countSources(cache.trials,groupMask);
  summaryRows(end+1,:)={groupOrder(g),domain,numel(find(groupMask)), ...
    numel(unique(session)),numel(unique(animal)),sourceCount, ...
    profiles{g}.adaptation,profiles{g}.noiseSDPC(1),profiles{g}.logLikelihood}; %#ok<AGROW>
  for r=1:height(P)
    profileRows(end+1,:)={groupOrder(g),domain,P.adaptation(r), ...
      P.logLikelihood(r),P.deltaLogLFromBest(r),P.noiseSDPC(r,1), ...
      P.evalCount(r,1)}; %#ok<AGROW>
  end
end
groupSummary=cell2table(summaryRows,'VariableNames', ...
  {'analysisGroup','domain','nTrials','nSessions','nAnimals','nSourceFiles', ...
   'bestAdaptation','noiseSDPC','logLikelihood'});
groupProfile=cell2table(profileRows,'VariableNames', ...
  {'analysisGroup','domain','adaptation','logLikelihood', ...
   'deltaLogLFromBest','noiseSDPC','evalCount'});
groupSummary.analysisGroup=string(groupSummary.analysisGroup);
groupSummary.domain=string(groupSummary.domain);
groupProfile.analysisGroup=string(groupProfile.analysisGroup);
groupProfile.domain=string(groupProfile.domain);

comparison=string([ ...
  "IDQ fixed versus randomized";"IDQ fixed versus randomized"; ...
  "Randomized IDQ versus IDS";"Randomized IDQ versus IDS"; ...
  "Randomized IDQ versus IDS"]);
model=string([ ...
  "pooled adaptation and noise";"separate adaptation and noise"; ...
  "pooled adaptation and noise";"shared adaptation, separate noise"; ...
  "separate adaptation and noise"]);
nParameters=[2;4;2;3;4];
logLikelihood=nan(5,1);
deltaLogLFromPooled=nan(5,1);
estimable=false(5,1);
note=strings(5,1);
idqPooled=[];
acrossSharedA=[];
acrossPooled=[];

hasIDQRandom=any(groupSummary.analysisGroup=="IDQ-randomized");
hasIDQFixed=any(groupSummary.analysisGroup=="IDQ-fixed");
if hasIDQRandom && hasIDQFixed
  idqMask=primaryMask & startsWith(trialGroup,"IDQ-");
  idqPooled=maxpoolfit.profileAdaptation(cache,cfg,'TrialMask',idqMask, ...
    'NDraws',p.Results.NDraws,'AdaptationGrid',p.Results.AdaptationGrid);
  idqSeparate=sum(groupSummary.logLikelihood(ismember( ...
    groupSummary.analysisGroup,["IDQ-randomized" "IDQ-fixed"])));
  logLikelihood(1:2)=[idqPooled.logLikelihood;idqSeparate];
  deltaLogLFromPooled(1:2)=[0;idqSeparate-idqPooled.logLikelihood];
  estimable(1:2)=true;
else
  missingIDQ="IDQ-randomized";
  if hasIDQRandom, missingIDQ="IDQ-fixed"; end
  note(1:2)="Not estimable: no primary trials for "+missingIDQ+".";
end

hasIDS=any(groupSummary.analysisGroup=="IDS");
if hasIDQRandom && hasIDS
  acrossMask=primaryMask & ismember(trialGroup,["IDQ-randomized" "IDS"]);
  acrossSharedA=maxpoolfit.profileAdaptation(cache,cfg,'TrialMask',acrossMask, ...
    'NDraws',p.Results.NDraws,'AdaptationGrid',p.Results.AdaptationGrid);
  pooledCache=cache;
  pooledCache.trials.domain(acrossMask)="IDQ_randomized_IDS";
  acrossPooled=maxpoolfit.profileAdaptation(pooledCache,cfg, ...
    'TrialMask',acrossMask,'NDraws',p.Results.NDraws, ...
    'AdaptationGrid',p.Results.AdaptationGrid);
  acrossSeparate=sum(groupSummary.logLikelihood(ismember( ...
    groupSummary.analysisGroup,["IDQ-randomized" "IDS"])));
  logLikelihood(3:5)=[acrossPooled.logLikelihood; ...
    acrossSharedA.logLikelihood;acrossSeparate];
  deltaLogLFromPooled(3:5)=[0; ...
    acrossSharedA.logLikelihood-acrossPooled.logLikelihood; ...
    acrossSeparate-acrossPooled.logLikelihood];
  estimable(3:5)=true;
else
  note(3:5)="Not estimable: randomized IDQ and IDS are both required.";
end
modelComparison=table(comparison,model,nParameters,logLikelihood, ...
  deltaLogLFromPooled,estimable,note);

result=struct();
result.expectedGroupOrder=expectedGroupOrder;
result.groupOrder=groupOrder;
result.groupProfiles=profiles;
result.groupSummary=groupSummary;
result.groupProfile=groupProfile;
result.modelComparison=modelComparison;
result.idqPooledProfile=idqPooled;
result.acrossSharedAdaptationProfile=acrossSharedA;
result.acrossFullyPooledProfile=acrossPooled;
result.trialGroup=trialGroup;
result.trialMask=primaryMask;
result.nNoiseDraws=p.Results.NDraws;
result.adaptationGrid=unique(double(p.Results.AdaptationGrid(:)),'stable');
result.config=cfg;
result.createdAt=datetime('now');
end

function n=countSources(T,mask)
if ismember('sourceFile',T.Properties.VariableNames)
  n=numel(unique(string(T.sourceFile(mask))));
else
  n=numel(unique(string(T.sessionID(mask))));
end
end
