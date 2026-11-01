function [pRDK1,details] = choiceProbabilityFromResponses(B,A,P,adaptation,noiseSD,cfg,z1,z2,returnWinners)
%CHOICEPROBABILITYFROMRESPONSES Fast likelihood core for precomputed terms.
if nargin < 8, returnWinners=false; end
assert(isequal(size(B),size(A),size(P)) && size(B,2)==2, ...
  'Response arrays must be nDirections x 2 x nTrials.');
assert(size(B,1)==size(z1,1) && isequal(size(z1),size(z2)), ...
  'Noise draws do not match response direction dimension.');
n=size(B,3); nDraws=size(z1,2); pRDK1=nan(n,1);
lambda=0;
if isfield(cfg,'independentNoiseFraction')
  lambda=double(cfg.independentNoiseFraction);
end
assert(isscalar(lambda)&&isfinite(lambda)&&lambda>=0&&lambda<=1, ...
  'cfg.independentNoiseFraction must be between zero and one.');
if lambda>0
  zi1=maxpoolfit.independentNoiseDraws(cfg,nDraws,2000003);
  zi2=maxpoolfit.independentNoiseDraws(cfg,nDraws,3000017);
  [noiseScale,~,~]=maxpoolfit.responseNoiseScale(B,cfg);
else
  zi1=[]; zi2=[]; noiseScale=[];
end
if returnWinners
  details.winnerMechanismIndexRDK1=zeros(n,nDraws,'uint16');
  details.winnerMechanismIndexRDK2=zeros(n,nDraws,'uint16');
else
  details=struct();
end
for first=1:cfg.batchSize:n
  q=first:min(n,first+cfg.batchSize-1);
  D=maxpoolfit.decisionResponse(B(:,:,q),A(:,:,q),P(:,:,q),adaptation,cfg);
  for b=1:numel(q)
    if lambda==0
      n1=noiseSD.*z1; n2=noiseSD.*z2;
    else
      n1=noiseSD.*(sqrt(1-lambda).*z1+ ...
        sqrt(lambda).*noiseScale(:,1,q(b)).*zi1);
      n2=noiseSD.*(sqrt(1-lambda).*z2+ ...
        sqrt(lambda).*noiseScale(:,2,q(b)).*zi2);
    end
    [m1,w1]=max(D(:,1,b)+n1,[],1);
    [m2,w2]=max(D(:,2,b)+n2,[],1);
    if cfg.useAllCrossPairs
      pRDK1(q(b))=crossPairProbability(m1,m2);
    else
      pRDK1(q(b))=mean(m1>m2)+0.5*mean(m1==m2);
    end
    if returnWinners
      details.winnerMechanismIndexRDK1(q(b),:)=uint16(w1);
      details.winnerMechanismIndexRDK2(q(b),:)=uint16(w2);
    end
  end
end
pRDK1=min(max(pRDK1,cfg.probabilityFloor),1-cfg.probabilityFloor);
end

function p=crossPairProbability(x,y)
% Exact Mann-Whitney comparison of all x-y pairs in O(n log n).
x=x(:); y=y(:); nX=numel(x); nY=numel(y);
[v,order]=sort([x;y]); isX=order<=nX;
score=0; nYLess=0; first=1;
while first<=numel(v)
  last=first;
  while last<numel(v) && v(last+1)==v(first), last=last+1; end
  group=isX(first:last);
  nXGroup=sum(group); nYGroup=numel(group)-nXGroup;
  score=score+nXGroup*(nYLess+0.5*nYGroup);
  nYLess=nYLess+nYGroup;
  first=last+1;
end
p=score/(nX*nY);
end
