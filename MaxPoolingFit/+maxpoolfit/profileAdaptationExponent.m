function result = profileAdaptationExponent(cache,cfg,varargin)
%PROFILEADAPTATIONEXPONENT Profile power weighting of prestep adaptation.
%
% The shared adaptation coefficient is held fixed while q is profiled in
%   Rstep - a*Ppre.*(Ppre/P0).^q
% and one internal-noise SD is reoptimized for each domain. q=0 is exactly
% the current linear prestep-adaptation model.

if ~isfield(cfg,'adaptationExponent'), cfg.adaptationExponent=0; end
if ~isfield(cfg,'adaptationReferencePC'), cfg.adaptationReferencePC=40; end
p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('ExponentGrid',0:0.1:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('Adaptation',cfg.adaptation, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
p.addParameter('ReferencePC',cfg.adaptationReferencePC, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(string(cache.schemaVersion)==string(cfg.cacheSchemaVersion), ...
  'Cache schema does not match the model configuration.');
if isfield(cache,'kernelStatus')
  assert(string(cache.kernelStatus)=="current", ...
    'Cached temporal predictors are stale; recompute kernels.');
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
exponentGrid=unique(double(p.Results.ExponentGrid(:)),'stable');
adaptation=double(p.Results.Adaptation);
referencePC=double(p.Results.ReferencePC);
domains=unique(string(cache.trials.domain(mask)),'stable');
assert(~isempty(domains),'No trials selected.');

z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
precomputed=cell(numel(domains),1);
for id=1:numel(domains)
  idx=find(mask & string(cache.trials.domain)==domains(id));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  precomputed{id}=struct('B',B,'A',A,'Rpre',Rpre, ...
    'choice1',double(cache.trials.choiceRDK(idx)==1));
end

nQ=numel(exponentGrid); nD=numel(domains);
noiseSDPC=nan(nQ,nD); domainLogLikelihood=nan(nQ,nD);
evalCount=zeros(nQ,nD);
for iq=1:nQ
  q=exponentGrid(iq);
  cfgQ=maxpoolfit.configureAdaptationExponent(cfg,q,referencePC);
  if cfg.verbose
    fprintf('Adaptation exponent %.3g (%d/%d)\n',q,iq,nQ);
  end
  for id=1:nD
    P=precomputed{id};
    objective=@(logSigma)nllAt(exp(logSigma),P,cfgQ);
    [x,fval,~,out]=fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
      log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
    noiseSDPC(iq,id)=exp(x);
    domainLogLikelihood(iq,id)=-fval;
    evalCount(iq,id)=out.funcCount;
  end
end
logLikelihood=sum(domainLogLikelihood,2);
[bestLogLikelihood,bestRow]=max(logLikelihood);
linearRow=find(abs(exponentGrid)<=1e-12,1);
if isempty(linearRow), deltaLogLFromLinear=nan(nQ,1);
else, deltaLogLFromLinear=logLikelihood-logLikelihood(linearRow); end
profile=table(exponentGrid,logLikelihood, ...
  logLikelihood-bestLogLikelihood,deltaLogLFromLinear,noiseSDPC, ...
  domainLogLikelihood,evalCount, ...
  'VariableNames',{'adaptationExponent','logLikelihood', ...
  'deltaLogLFromBest','deltaLogLFromLinear','noiseSDPC', ...
  'domainLogLikelihood','evalCount'});
bestCfg=maxpoolfit.configureAdaptationExponent(cfg,exponentGrid(bestRow),referencePC);
linearCfg=maxpoolfit.configureAdaptationExponent(cfg,0,referencePC);
result=struct('modelVersion',bestCfg.modelVersion, ...
  'baseModelVersion',linearCfg.modelVersion, ...
  'cacheSchemaVersion',cache.schemaVersion,'domains',domains, ...
  'profile',profile,'bestRow',bestRow, ...
  'adaptationExponent',exponentGrid(bestRow), ...
  'adaptation',adaptation,'adaptationReferencePC',referencePC, ...
  'logLikelihood',bestLogLikelihood, ...
  'noiseSDPC',noiseSDPC(bestRow,:).','nTrials',sum(mask), ...
  'trialMask',mask,'nNoiseDraws',p.Results.NDraws, ...
  'config',bestCfg,'createdAt',datetime('now'));

  function nll=nllAt(sigma,P,cfgQ)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      P.B,P.A,P.Rpre,adaptation,sigma,cfgQ,z1,z2,false);
    pChosen=P.choice1.*p1+(1-P.choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end
end
