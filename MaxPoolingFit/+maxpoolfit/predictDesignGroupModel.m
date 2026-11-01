function prediction = predictDesignGroupModel(cache,fit,cfg,varargin)
%PREDICTDESIGNGROUPMODEL Predict with group-specific adaptation and noise.
p=inputParser;
p.addParameter('TrialMask',fit.trialMask,@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',fit.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('ReturnWinners',false,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
if islogical(p.Results.TrialMask), idx=find(p.Results.TrialMask);
else, idx=p.Results.TrialMask(:); end
pRDK1=nan(height(cache.trials),1);
winnerDetails=cell(height(fit.groupFit),1);
for g=1:height(fit.groupFit)
  q=idx(fit.trialGroup(idx)==fit.groupFit.analysisGroup(g));
  [p1,d]=maxpoolfit.choiceProbability(cache,q,fit.groupFit.adaptation(g), ...
    fit.groupFit.noiseSDPC(g),cfg,'NDraws',p.Results.NDraws, ...
    'ReturnWinners',p.Results.ReturnWinners);
  pRDK1(q)=p1;
  winnerDetails{g}=struct('trialIdx',q,'details',d);
end
assert(all(isfinite(pRDK1(idx))),'One or more selected trials were not predicted.');
prediction=struct('pRDK1',pRDK1,'trialIdx',idx, ...
  'adaptation',NaN,'adaptationByGroup',fit.adaptationByGroup, ...
  'trialGroup',fit.trialGroup,'groupFit',fit.groupFit, ...
  'nNoiseDraws',p.Results.NDraws,'winnerDetails',{winnerDetails}, ...
  'createdAt',datetime('now'));
end
