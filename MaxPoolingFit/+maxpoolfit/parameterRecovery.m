function recovery = parameterRecovery(baseCache, truth, cfg, varargin)
%PARAMETERRECOVERY Simulate and refit repeated datasets.
p=inputParser;
p.addParameter('NReplicates',10,@(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('Seed',2200,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
n=p.Results.NReplicates; ll=nan(n,1);
domains=unique(string(baseCache.trials.domain),'stable'); sigmaHat=nan(n,numel(domains));
fits=cell(n,1);
for r=1:n
  sim=maxpoolfit.simulateChoices(baseCache,cfg.adaptation,truth.noiseSDPC,cfg, ...
    'Seed',p.Results.Seed+r);
  fits{r}=maxpoolfit.fitModel(sim,cfg);
  sigmaHat(r,:)=fits{r}.noiseSDPC(:)'; ll(r)=fits{r}.logLikelihood;
end
recovery=struct('truth',truth,'domains',domains, ...
  'noiseSDPChat',sigmaHat,'logLikelihood',ll,'fits',{fits}, ...
  'nReplicates',n,'createdAt',datetime('now'));
end
