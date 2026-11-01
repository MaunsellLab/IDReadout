function prediction = predict(cache, fit, cfg, varargin)
%PREDICT Fitted probabilities and optional winner-direction diagnostics.
p=inputParser;
p.addParameter('TrialMask',true(height(cache.trials),1),@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.finalNoiseDraws,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('ReturnWinners',false,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
if ~isfield(cfg,'adaptationExponent'), cfg.adaptationExponent=0; end
if ~isfield(cfg,'adaptationGainStrength'), cfg.adaptationGainStrength=0; end
if ~isfield(cfg,'adaptationReferencePC'), cfg.adaptationReferencePC=40; end
if islogical(p.Results.TrialMask), idx=find(p.Results.TrialMask); else, idx=p.Results.TrialMask(:); end
pRDK1=nan(height(cache.trials),1); detailsByDomain=cell(numel(fit.domains),1);
for k=1:numel(fit.domains)
  q=idx(string(cache.trials.domain(idx))==fit.domains(k));
  [p1,d]=maxpoolfit.choiceProbability(cache,q,fit.adaptation,fit.noiseSDPC(k),cfg, ...
    'NDraws',p.Results.NDraws,'ReturnWinners',p.Results.ReturnWinners);
  pRDK1(q)=p1; detailsByDomain{k}=struct('trialIdx',q,'details',d);
end
prediction=struct('pRDK1',pRDK1,'trialIdx',idx,'adaptation',fit.adaptation, ...
  'adaptationExponent',cfg.adaptationExponent, ...
  'adaptationGainStrength',cfg.adaptationGainStrength, ...
  'adaptationReferencePC',cfg.adaptationReferencePC, ...
  'domains',fit.domains,'noiseSDPC',fit.noiseSDPC,'nNoiseDraws',p.Results.NDraws, ...
  'winnerDetails',{detailsByDomain},'createdAt',datetime('now'));
end
