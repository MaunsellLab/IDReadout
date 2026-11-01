function result = profileAdaptation(cache,cfg,varargin)
%PROFILEADAPTATION Profile corrected adaptation with domain-specific noise.
%
% At each tested value of a, one
% internal-noise SD is reoptimized for each experimental domain.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('AdaptationGrid',0:0.1:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(string(cache.schemaVersion)==string(cfg.cacheSchemaVersion), ...
  'Cache schema does not match the model configuration.');
if isfield(cache,'kernelStatus')
  assert(string(cache.kernelStatus)=="current", ...
    'Cached temporal predictors are stale; rerun cache synchronization with RecomputeKernels=true.');
end
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
adaptationGrid=unique(double(p.Results.AdaptationGrid(:)),'stable');
domains=unique(string(cache.trials.domain(mask)),'stable');
assert(~isempty(domains),'No trials selected.');

z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
precomputed=cell(numel(domains),1);
for id=1:numel(domains)
  idx=find(mask & string(cache.trials.domain)==domains(id));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  precomputed{id}=struct('idx',idx,'B',B,'A',A,'Rpre',Rpre, ...
    'choice1',double(cache.trials.choiceRDK(idx)==1));
end

nA=numel(adaptationGrid); nD=numel(domains);
noiseSDPC=nan(nA,nD);
domainLogLikelihood=nan(nA,nD);
evalCount=zeros(nA,nD);
for ia=1:nA
  a=adaptationGrid(ia);
  if cfg.verbose
    fprintf('Adaptation %.3g (%d/%d)\n',a,ia,nA);
  end
  for id=1:nD
    P=precomputed{id};
    objective=@(logSigma)nllAt(exp(logSigma),a,P);
    [x,fval,~,out]=fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
      log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
    noiseSDPC(ia,id)=exp(x);
    domainLogLikelihood(ia,id)=-fval;
    evalCount(ia,id)=out.funcCount;
  end
end
logLikelihood=sum(domainLogLikelihood,2);
[bestLogLikelihood,bestRow]=max(logLikelihood);
zeroRow=find(abs(adaptationGrid)<=1e-12,1);
if isempty(zeroRow)
  deltaLogLFromZero=nan(nA,1);
else
  deltaLogLFromZero=logLikelihood-logLikelihood(zeroRow);
end
profile=table(adaptationGrid,logLikelihood, ...
  logLikelihood-bestLogLikelihood,deltaLogLFromZero,noiseSDPC, ...
  domainLogLikelihood,evalCount, ...
  'VariableNames',{'adaptation','logLikelihood','deltaLogLFromBest', ...
  'deltaLogLFromZero','noiseSDPC','domainLogLikelihood','evalCount'});
result=struct();
result.modelVersion=cfg.modelVersion;
result.cacheSchemaVersion=cache.schemaVersion;
result.domains=domains;
result.profile=profile;
result.bestRow=bestRow;
result.adaptation=adaptationGrid(bestRow);
result.logLikelihood=bestLogLikelihood;
result.noiseSDPC=noiseSDPC(bestRow,:).';
result.nTrials=sum(mask);
result.trialMask=mask;
result.nNoiseDraws=p.Results.NDraws;
result.config=cfg;
result.createdAt=datetime('now');

  function nll=nllAt(sigma,a,P)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      P.B,P.A,P.Rpre,a,sigma,cfg,z1,z2,false);
    pChosen=P.choice1.*p1+(1-P.choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end
end
