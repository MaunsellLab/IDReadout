function [pRDK1,audit] = predict(R,p,noiseSD,Z)
%PREDICT Compare signed post-pooling changes in the two RDKs.
% Gaussian fluctuations are added to absolute rates before pooling. Negative
% absolute rates are floored at zero (a firing-rate bound). Neither unit
% pre/post changes nor pooled pre/post differences are rectified.
assert(noiseSD>=0 && isfinite(noiseSD),'Invalid internal-noise SD.');
assert(size(R,2)==2 && size(R,3)==2,'Expected two RDKs and two epochs.');
n=size(R,4); nDraws=size(Z,2); pRDK1=zeros(n,1);
clipped=0; evaluated=0;
for k=1:n
  delta=zeros(nDraws,2);
  for j=1:2
    S=zeros(nDraws,2);
    for e=1:2
      x=R(:,j,e,k)+noiseSD*Z(:,:,j,e);
      clipped=clipped+nnz(x<0); evaluated=evaluated+numel(x);
      S(:,e)=mtpnorm.pool(max(x,0),p)';
    end
    delta(:,j)=S(:,2)-S(:,1);
  end
  % Integrate all independent RDK draw pairs, including half-weighted ties.
  x=delta(:,1); y=sort(delta(:,2));
  % Sorting/searching below handles repeated y values and exact ties.
  [v,order]=sort([x;y]); isX=order<=nDraws;
  score=0; nYLess=0; first=1;
  while first<=numel(v)
    last=first;
    while last<numel(v) && v(last+1)==v(first), last=last+1; end
    g=isX(first:last); nx=sum(g); ny=numel(g)-nx;
    score=score+nx*(nYLess+.5*ny); nYLess=nYLess+ny; first=last+1;
  end
  pRDK1(k)=score/nDraws^2;
end
audit=struct('negativeRateFraction',clipped/max(evaluated,1), ...
  'nDrawsPerRDK',nDraws,'nIndependentRDKPairs',nDraws^2, ...
  'pooledDeltaRectified',false,'unitDeltaRectified',false);
end
