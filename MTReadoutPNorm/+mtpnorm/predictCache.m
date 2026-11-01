function [prob,audit,responseAudit] = predictCache(cache,cfg,p,sd,Z,idx)
%PREDICTCACHE Bound response-array memory independently of session count.
if nargin<6, idx=(1:height(cache.trials))'; end
if islogical(idx), idx=find(idx); end
idx=idx(:); n=numel(idx); prob=zeros(n,1); clipped=0; evaluated=0;
assert(n>0,'No trials to predict.');
maxBudget=0; minCoh=Inf; responseAudit=struct();
for first=1:cfg.trialBatchSize:n
  q=first:min(n,first+cfg.trialBatchSize-1);
  [R,A]=mtpnorm.fullResponses(cache,cfg,idx(q));
  [prob(q),a]=mtpnorm.predict(R,p,sd,Z);
  count=numel(R)*size(Z,2);
  clipped=clipped+a.negativeRateFraction*count; evaluated=evaluated+count;
  maxBudget=max(maxBudget,A.maxCoherentBudgetPC);
  minCoh=min(minCoh,A.minStreamCoherencePC); responseAudit=A;
end
audit=a; audit.negativeRateFraction=clipped/max(evaluated,1);
responseAudit.maxCoherentBudgetPC=maxBudget;
responseAudit.minStreamCoherencePC=minCoh;
end
