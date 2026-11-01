function [baseResponse, adaptationTemplate, preResponse] = deterministicResponse(cache, idx, cfg)
%DETERMINISTICRESPONSE Construct step, adaptation, and prestep responses.
theta = cfg.directionGridDeg(:);
n = numel(idx); m = numel(theta);
baseResponse = zeros(m,2,n,'double');
adaptationTemplate = zeros(m,2,n,'double');
preResponse = zeros(m,2,n,'double');
for q = 1:n
  t = idx(q);
  for j = 1:2
    dirs = double(cache.streamDirDeg{t,j}(:));
    cohs = double(cache.stepCohPC{t,j}(:));
    R = maxpoolfit.coherenceTuning(theta,dirs,cfg.sigmaMTDeg,cfg) * cohs;
    preCohs = double(cache.preCohPC{t,j}(:));
    T = maxpoolfit.coherenceTuning(theta,dirs,cfg.sigmaAdaptDeg,cfg) * preCohs;
    P = maxpoolfit.coherenceTuning(theta,dirs,cfg.sigmaMTDeg,cfg) * preCohs;
    baseResponse(:,j,q) = R;
    adaptationTemplate(:,j,q) = T;
    preResponse(:,j,q) = P;
  end
end
end
