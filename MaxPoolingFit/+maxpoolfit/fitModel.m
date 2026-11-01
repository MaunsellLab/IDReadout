function fit = fitModel(cache, cfg, varargin)
%FITMODEL Fit domain-noise SDs at a specified shared adaptation value.
p = inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
if ~isfield(cfg,'adaptationExponent'), cfg.adaptationExponent=0; end
if ~isfield(cfg,'adaptationGainStrength'), cfg.adaptationGainStrength=0; end
if ~isfield(cfg,'adaptationReferencePC'), cfg.adaptationReferencePC=40; end
maxpoolfit.validateCache(cache);
assert(string(cache.schemaVersion)==string(cfg.cacheSchemaVersion), ...
  'Cache schema does not match the model configuration.');
if isfield(cache,'kernelStatus')
  assert(string(cache.kernelStatus)=="current", ...
    'Cached temporal predictors are stale; rerun cache synchronization with RecomputeKernels=true.');
end
if isempty(p.Results.TrialMask)
  mask = cfg.primarySelector(cache.trials);
else
  mask = false(height(cache.trials),1);
  if islogical(p.Results.TrialMask), mask=p.Results.TrialMask(:); else, mask(p.Results.TrialMask)=true; end
end
domains = unique(string(cache.trials.domain(mask)),'stable');
assert(~isempty(domains),'No trials selected.');
assert(isfield(cfg,'adaptation') && isnumeric(cfg.adaptation) && ...
  isscalar(cfg.adaptation) && isfinite(cfg.adaptation) && ...
  cfg.adaptation>=0 && cfg.adaptation<=1, ...
  'cfg.adaptation must be a finite scalar between zero and one.');

precomputed=cell(numel(domains),1);
z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
for id=1:numel(domains)
  idx=find(mask & string(cache.trials.domain)==domains(id));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  precomputed{id}=struct('idx',idx,'B',B,'A',A,'Rpre',Rpre, ...
    'choice1',double(cache.trials.choiceRDK(idx)==1));
end

sigmaHat = nan(numel(domains),1);
llDomain = nan(numel(domains),1);
evalCount = zeros(numel(domains),1);
nTrialsDomain = zeros(numel(domains),1);
for id = 1:numel(domains)
  if cfg.verbose, fprintf('Domain %s (%d/%d)\n',domains(id),id,numel(domains)); end
  P=precomputed{id};
  objective = @(logSigma) nllAt(exp(logSigma),P);
  [x,fval,~,out] = fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
    log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
  sigmaHat(id)=exp(x);
  llDomain(id)=-fval;
  evalCount(id)=out.funcCount;
  nTrialsDomain(id)=numel(P.idx);
end
fit = struct();
fit.modelVersion=cfg.modelVersion; fit.cacheSchemaVersion=cache.schemaVersion;
fit.adaptation=cfg.adaptation; fit.domains=domains;
fit.adaptationExponent=cfg.adaptationExponent;
fit.adaptationGainStrength=cfg.adaptationGainStrength;
fit.adaptationReferencePC=cfg.adaptationReferencePC;
fit.noiseSDPC=sigmaHat; fit.logLikelihood=sum(llDomain);
fit.domainFit=table(domains,nTrialsDomain,sigmaHat,llDomain,evalCount, ...
  'VariableNames',{'domain','nTrials','noiseSDPC','logLikelihood','evalCount'});
fit.nTrials=sum(mask); fit.trialMask=mask; fit.nNoiseDraws=p.Results.NDraws;
fit.config=cfg; fit.createdAt=datetime('now');

  function nll = nllAt(sigma,P)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      P.B,P.A,P.Rpre,cfg.adaptation,sigma,cfg,z1,z2,false);
    pChosen=P.choice1.*p1+(1-P.choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end
end
