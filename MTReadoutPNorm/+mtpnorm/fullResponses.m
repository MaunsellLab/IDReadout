function [R,audit] = fullResponses(cache,cfg,idx)
%FULLRESPONSES Expected full response to fixed-budget direction mixtures.
% R: direction x physical RDK x {pre,step} x trial. The remaining dots
% are assumed isotropically directed; exact realized random-dot angles are
% not present in the coherence traces. Integration PRECEDES pooling.
if nargin<3, idx=(1:height(cache.trials))'; end
if islogical(idx), idx=find(idx); end
theta=(0:cfg.directionSpacingDeg:360-cfg.directionSpacingDeg)';
assert(abs(numel(theta)*cfg.directionSpacingDeg-360)<1e-8, ...
  'Direction spacing must divide 360.');
assert(cfg.spontaneousResponse>=0,'Spontaneous response must be nonnegative.');
% Continuous circular-Gaussian mean; independent of the response-bank grid.
gbar=sqrt(2*pi)*cfg.sigmaMTDeg/360*erf(180/(sqrt(2)*cfg.sigmaMTDeg));
n=numel(idx); R=zeros(numel(theta),2,2,n);
largestBudget=0; smallestCoherence=Inf;
for k=1:n
  source=idx(k);
  dirs=cache.streamDirDeg{source}(:)';
  d=mod(theta-dirs+180,360)-180;
  G=exp(-.5*(d/cfg.sigmaMTDeg).^2);
  for j=1:2
    cohs={cache.preCohPC{source,j},cache.stepCohPC{source,j}};
    for e=1:2
      c=double(cohs{e}(:));
      assert(numel(c)==size(G,2),'Stream count mismatch.');
      assert(all(isfinite(c)) && min(c)>=-1e-4 && sum(c)<=100+1e-4, ...
        'Unphysical dot budget in trial %d, RDK %d, epoch %d.',source,j,e);
      smallestCoherence=min(smallestCoherence,min(c));
      c=max(c,0); largestBudget=max(largestBudget,sum(c));
      R(:,j,e,k)=cfg.spontaneousResponse+(100-sum(c))*gbar+G*c;
    end
  end
end
audit=struct('randomDirectionMean',gbar,'maxCoherentBudgetPC',largestBudget, ...
  'minStreamCoherencePC',smallestCoherence,'nDirections',numel(theta), ...
  'randomDotRemainder',"isotropic expected distribution");
end
