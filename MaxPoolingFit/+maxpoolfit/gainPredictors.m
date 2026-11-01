function P = gainPredictors(cache,trialIdx)
%GAINPREDICTORS Reconstruct changed/no-change drift and probe noise inputs.
%
% Predictors are kernel-weighted physical coherence values. Drift predictors
% have their known deterministic coherence removed. Probe pedestal constants
% are retained here and removed by within-session centering in gainComparison.
% For IDR, one predictor represents the single noise trace yoked across the
% displayed probe pair. For IDQ/IDS, the two independent probe predictors are
% summed so their constrained coefficient is the gain per probe stream.

maxpoolfit.validateCache(cache);
trialIdx=double(trialIdx(:));
assert(all(trialIdx>=1 & trialIdx<=height(cache.trials)), ...
  'Trial indices are outside the cache.');
T=cache.trials(trialIdx,:); n=numel(trialIdx);

changedDriftNoisePC=nan(n,1); unchangedDriftNoisePC=nan(n,1);
changedProbeNoisePC=nan(n,1); unchangedProbeNoisePC=nan(n,1);
probePredictorContract=strings(n,1);
for k=1:n
  c=double(T.changeRDK(k)); u=3-c;
  changed=double(cache.stepCohPC{trialIdx(k),c}(:));
  unchanged=double(cache.stepCohPC{trialIdx(k),u}(:));
  assert(numel(changed)>=2 && numel(unchanged)>=2, ...
    'Every gain-analysis trial must contain drift and probe streams.');
  changedDriftNoisePC(k)=changed(1)-double(T.stepEndpointCohPC(k));
  unchangedDriftNoisePC(k)=unchanged(1)-double(T.baseCohPC(k));
  if string(T.domain(k))=="IDR"
    changedProbeNoisePC(k)=changed(2);
    unchangedProbeNoisePC(k)=unchanged(2);
    probePredictorContract(k)="one yoked probe-noise trace";
  else
    changedProbeNoisePC(k)=sum(changed(2:end));
    unchangedProbeNoisePC(k)=sum(unchanged(2:end));
    probePredictorContract(k)="sum of two independent probe-noise traces";
  end
end

P=T(:,intersect({'domain','animal','sessionID','trialID','changeRDK', ...
  'stepMagnitudePC','probeOffsetDeg'},T.Properties.VariableNames,'stable'));
P.cacheTrialIdx=trialIdx;
P.changedDriftNoisePC=changedDriftNoisePC;
P.unchangedDriftNoisePC=unchangedDriftNoisePC;
P.changedProbeNoisePC=changedProbeNoisePC;
P.unchangedProbeNoisePC=unchangedProbeNoisePC;
P.probePredictorContract=probePredictorContract;
end
