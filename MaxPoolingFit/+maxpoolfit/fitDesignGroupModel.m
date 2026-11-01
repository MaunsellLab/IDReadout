function fit = fitDesignGroupModel(cache,cfg,trialGroup,adaptationByGroup,varargin)
%FITDESIGNGROUPMODEL Fit one noise SD per design group at fixed adaptations.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
trialGroup=string(trialGroup(:));
assert(numel(trialGroup)==height(cache.trials), ...
  'trialGroup must contain one label per cached trial.');
assert(istable(adaptationByGroup)&&all(ismember( ...
  {'analysisGroup','adaptation'},adaptationByGroup.Properties.VariableNames)), ...
  'adaptationByGroup must contain analysisGroup and adaptation.');
adaptationByGroup.analysisGroup=string(adaptationByGroup.analysisGroup);
assert(numel(unique(adaptationByGroup.analysisGroup))==height(adaptationByGroup), ...
  'adaptationByGroup labels must be unique.');
assert(all(isfinite(adaptationByGroup.adaptation))&& ...
  all(adaptationByGroup.adaptation>=0 & adaptationByGroup.adaptation<=1), ...
  'All design-group adaptation values must be between zero and one.');

if isempty(p.Results.TrialMask)
  mask=cfg.primarySelector(cache.trials);
else
  mask=false(height(cache.trials),1);
  if islogical(p.Results.TrialMask)
    assert(numel(p.Results.TrialMask)==height(cache.trials), ...
      'Logical TrialMask must have one element per cached trial.');
    mask=p.Results.TrialMask(:);
  else
    mask(p.Results.TrialMask)=true;
  end
end
groups=adaptationByGroup.analysisGroup(ismember( ...
  adaptationByGroup.analysisGroup,unique(trialGroup(mask),'stable')));
assert(~isempty(groups),'No design-group trials were selected.');
assert(all(ismember(unique(trialGroup(mask)),groups)), ...
  'One or more selected trial groups lack a specified adaptation value.');

z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
nGroups=numel(groups);
nTrials=zeros(nGroups,1); nSessions=zeros(nGroups,1);
adaptation=zeros(nGroups,1); noiseSDPC=nan(nGroups,1);
logLikelihood=nan(nGroups,1); evalCount=zeros(nGroups,1);
domain=strings(nGroups,1);
for g=1:nGroups
  idx=find(mask & trialGroup==groups(g));
  nTrials(g)=numel(idx);
  nSessions(g)=numel(unique(string(cache.trials.sessionID(idx))));
  d=unique(string(cache.trials.domain(idx)),'stable');
  assert(numel(d)==1,'Each analysis group must contain one source domain.');
  domain(g)=d;
  adaptation(g)=adaptationByGroup.adaptation( ...
    adaptationByGroup.analysisGroup==groups(g));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  choice1=double(cache.trials.choiceRDK(idx)==1);
  objective=@(logSigma)nllAt(exp(logSigma),adaptation(g),B,A,Rpre,choice1);
  if cfg.verbose
    fprintf('Design group %s (%d/%d; a = %.3g)\n', ...
      groups(g),g,nGroups,adaptation(g));
  end
  [x,fval,~,out]=fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
    log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
  noiseSDPC(g)=exp(x); logLikelihood(g)=-fval;
  evalCount(g)=out.funcCount;
end

groupFit=table(groups,domain,nTrials,nSessions,adaptation,noiseSDPC, ...
  logLikelihood,evalCount,'VariableNames',{'analysisGroup','domain', ...
  'nTrials','nSessions','adaptation','noiseSDPC','logLikelihood','evalCount'});
trialAdaptation=nan(height(cache.trials),1);
for g=1:nGroups, trialAdaptation(trialGroup==groups(g))=adaptation(g); end
fit=struct();
fit.modelVersion=cfg.modelVersion;
fit.cacheSchemaVersion=cache.schemaVersion;
fit.adaptation=NaN;
fit.adaptationByGroup=groupFit(:,{'analysisGroup','adaptation'});
fit.groups=groups;
fit.groupFit=groupFit;
fit.noiseSDPC=noiseSDPC;
fit.logLikelihood=sum(logLikelihood);
fit.nTrials=sum(mask);
fit.trialMask=mask;
fit.trialGroup=trialGroup;
fit.trialAdaptation=trialAdaptation;
fit.nNoiseDraws=p.Results.NDraws;
fit.config=cfg;
fit.createdAt=datetime('now');

  function nll=nllAt(sigma,a,B,A,Rpre,choice1)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      B,A,Rpre,a,sigma,cfg,z1,z2,false);
    pChosen=choice1.*p1+(1-choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end
end
