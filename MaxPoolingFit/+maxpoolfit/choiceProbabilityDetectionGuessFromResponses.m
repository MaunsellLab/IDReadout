function [pRDK1,details] = choiceProbabilityDetectionGuessFromResponses( ...
    B,A,Rpre,adaptation,noiseSD,kappa,cfg,z1,z2,guessP_RDK1)
%CHOICEPROBABILITYDETECTIONGUESSFROMRESPONSES Detection race with guessing.
%
% Each RDK's noisy step maximum is centered on its deterministic no-change
% maximum. If neither centered maximum reaches kappa*noiseSD, the model
% guesses. Otherwise it retains the original absolute-max RDK comparison.
% kappa=-Inf is the exact original absolute-max comparison.

if nargin<11, guessP_RDK1=0.5; end
assert(isequal(size(B),size(A),size(Rpre)) && size(B,2)==2, ...
  'Step, prestep, and adaptation arrays must match and contain two RDKs.');
assert(isfinite(noiseSD) && noiseSD>0,'noiseSD must be positive and finite.');
assert(isscalar(kappa) && (isfinite(kappa)||isinf(kappa)), ...
  'kappa must be a scalar.');
assert(kappa>=0 || (isinf(kappa)&&kappa<0), ...
  'kappa must be nonnegative or -Inf for the no-guess model.');
assert(isfinite(guessP_RDK1) && guessP_RDK1>=0 && guessP_RDK1<=1, ...
  'guessP_RDK1 must lie between zero and one.');

n=size(B,3); pRDK1=nan(n,1); guessProbability=zeros(n,1);
D=maxpoolfit.decisionResponse(B,A,Rpre,adaptation,cfg);
D0=maxpoolfit.decisionResponse(Rpre,A,Rpre,adaptation,cfg);
for q=1:n
  m1=max(D(:,1,q)+noiseSD.*z1,[],1);
  m2=max(D(:,2,q)+noiseSD.*z2,[],1);
  if isinf(kappa) && kappa<0
    pRDK1(q)=crossPairProbability(m1,m2);
    continue;
  end
  baseline1=max(D0(:,1,q)); baseline2=max(D0(:,2,q));
  c1=m1-baseline1; c2=m2-baseline2;
  threshold=kappa.*noiseSD;
  below1=c1<threshold; below2=c2<threshold;
  guessProbability(q)=mean(below1).*mean(below2);
  denominator=numel(m1).*numel(m2);
  originalWin=crossPairScore(m1,m2)./denominator;
  belowWin=crossPairScore(m1(below1),m2(below2))./denominator;
  pRDK1(q)=originalWin-belowWin+ ...
    guessP_RDK1.*guessProbability(q);
end
pRDK1=min(max(pRDK1,cfg.probabilityFloor),1-cfg.probabilityFloor);
details=struct('guessProbability',guessProbability,'kappa',kappa, ...
  'guessP_RDK1',guessP_RDK1,'meanGuessProbability',mean(guessProbability));
end

function p=crossPairProbability(x,y)
p=crossPairScore(x,y)./(numel(x).*numel(y));
end

function score=crossPairScore(x,y)
% Mann-Whitney score: number of x>y pairs plus half the tied pairs.
x=x(:); y=y(:); nX=numel(x); nY=numel(y);
if nX==0 || nY==0, score=0; return; end
[v,order]=sort([x;y]); isX=order<=nX;
score=0; nYLess=0; first=1;
while first<=numel(v)
  last=first;
  while last<numel(v) && v(last+1)==v(first), last=last+1; end
  group=isX(first:last);
  nXGroup=sum(group); nYGroup=numel(group)-nXGroup;
  score=score+nXGroup.*(nYLess+0.5.*nYGroup);
  nYLess=nYLess+nYGroup;
  first=last+1;
end
end
