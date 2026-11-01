function cache = simulateChoices(cache, adaptation, noiseByDomain, cfg, varargin)
%SIMULATECHOICES Draw choices from the model for parameter recovery.
p=inputParser; p.addParameter('Seed',7717,@(x)isnumeric(x)&&isscalar(x)); p.parse(varargin{:});
domains=unique(string(cache.trials.domain),'stable');
assert(numel(noiseByDomain)==numel(domains),'One noise SD is required per domain.');
p1=nan(height(cache.trials),1);
for k=1:numel(domains)
  idx=find(string(cache.trials.domain)==domains(k));
  p1(idx)=maxpoolfit.choiceProbability(cache,idx,adaptation,noiseByDomain(k),cfg);
end
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(p.Results.Seed,'twister');
cache.trials.choiceRDK=uint8(1+(rand(height(cache.trials),1)>=p1));
cache.trials.correct=cache.trials.choiceRDK==cache.trials.changeRDK;
cache.simulationTruth=struct('adaptation',adaptation,'domains',domains, ...
  'noiseSDPC',noiseByDomain,'seed',p.Results.Seed);
end
