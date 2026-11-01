function result = profileGuessCutoff(cache,cfg,varargin)
%PROFILEGUESSCUTOFF Profile a shared detection cutoff with domain noise.
%
% kappa=-Inf is the exact no-guess model. Finite kappa defines a detection
% threshold of kappa*noiseSD within every domain. Adaptation is held fixed.

p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('KappaGrid',[-Inf 0:0.25:3], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&& ...
  all((isfinite(x)&x>=0)|(isinf(x)&x<0)));
p.addParameter('GuessP_RDK1',0.5, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(string(cache.schemaVersion)==string(cfg.cacheSchemaVersion), ...
  'Cache schema does not match the model configuration.');
if isfield(cfg,'decisionRule')
  assert(string(cfg.decisionRule)=="absolute-step-max", ...
    'Guess-cutoff profiling must start from the absolute-step-max model.');
end
if isempty(p.Results.TrialMask)
  mask=cfg.primarySelector(cache.trials);
else
  mask=false(height(cache.trials),1);
  if islogical(p.Results.TrialMask)
    assert(numel(p.Results.TrialMask)==height(cache.trials), ...
      'Logical TrialMask must have one element per cached trial.');
    mask=p.Results.TrialMask(:);
  else
    mask(p.Results.TrialMask)=true;
  end
end
kappaGrid=unique(double(p.Results.KappaGrid(:)),'stable');
assert(any(isinf(kappaGrid)&kappaGrid<0), ...
  'KappaGrid must include -Inf as the no-guess reference.');
domains=unique(string(cache.trials.domain(mask)),'stable');
assert(~isempty(domains),'No trials selected.');

z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
precomputed=cell(numel(domains),1);
for id=1:numel(domains)
  idx=find(mask & string(cache.trials.domain)==domains(id));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  precomputed{id}=struct('idx',idx,'B',B,'A',A,'Rpre',Rpre, ...
    'choice1',double(cache.trials.choiceRDK(idx)==1), ...
    'change1',double(cache.trials.changeRDK(idx)==1));
end

nK=numel(kappaGrid); nD=numel(domains);
noiseSDPC=nan(nK,nD); domainLogLikelihood=nan(nK,nD);
guessFraction=nan(nK,nD); predictedAccuracy=nan(nK,nD);
evalCount=zeros(nK,nD);
for ik=1:nK
  kappa=kappaGrid(ik);
  if cfg.verbose
    if isinf(kappa), label='no guess'; else, label=sprintf('kappa %.3g',kappa); end
    fprintf('Guess cutoff %s (%d/%d)\n',label,ik,nK);
  end
  for id=1:nD
    P=precomputed{id}; objective=@(logSigma)nllAt(exp(logSigma),kappa,P);
    [x,fval,~,out]=fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
      log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
    sigma=exp(x); noiseSDPC(ik,id)=sigma;
    domainLogLikelihood(ik,id)=-fval; evalCount(ik,id)=out.funcCount;
    [p1,details]=probabilityAt(sigma,kappa,P);
    guessFraction(ik,id)=mean(details.guessProbability);
    predictedAccuracy(ik,id)=mean(P.change1.*p1+(1-P.change1).*(1-p1));
  end
end
logLikelihood=sum(domainLogLikelihood,2);
[bestLogLikelihood,bestRow]=max(logLikelihood);
referenceRow=find(isinf(kappaGrid)&kappaGrid<0,1);
deltaLogLFromNoGuess=logLikelihood-logLikelihood(referenceRow);
profile=table(kappaGrid,logLikelihood,logLikelihood-bestLogLikelihood, ...
  deltaLogLFromNoGuess,noiseSDPC,domainLogLikelihood,guessFraction, ...
  predictedAccuracy,evalCount, ...
  'VariableNames',{'kappa','logLikelihood','deltaLogLFromBest', ...
  'deltaLogLFromNoGuess','noiseSDPC','domainLogLikelihood', ...
  'guessFraction','predictedAccuracy','evalCount'});

bestP_RDK1=nan(height(cache.trials),1);
bestGuessProbability=nan(height(cache.trials),1);
for id=1:nD
  P=precomputed{id};
  [p1,details]=probabilityAt(noiseSDPC(bestRow,id),kappaGrid(bestRow),P);
  bestP_RDK1(P.idx)=p1;
  bestGuessProbability(P.idx)=details.guessProbability;
end
conditionSummary=conditionDiagnostics(cache,mask,bestP_RDK1, ...
  bestGuessProbability);

result=struct('modelVersion',"DenseHardMax-v7-detection-guess", ...
  'baseModelVersion',string(cfg.modelVersion), ...
  'cacheSchemaVersion',cache.schemaVersion,'adaptation',cfg.adaptation, ...
  'domains',domains,'profile',profile,'bestRow',bestRow, ...
  'kappa',kappaGrid(bestRow),'logLikelihood',bestLogLikelihood, ...
  'noiseSDPC',noiseSDPC(bestRow,:).','nTrials',sum(mask), ...
  'trialMask',mask,'nNoiseDraws',p.Results.NDraws, ...
  'guessP_RDK1',p.Results.GuessP_RDK1, ...
  'conditionSummary',conditionSummary, ...
  'prediction',struct('pRDK1',bestP_RDK1,'trialIdx',find(mask), ...
    'guessProbability',bestGuessProbability), ...
  'config',cfg,'createdAt',datetime('now'));

  function nll=nllAt(sigma,kappa,P)
    p1=probabilityAt(sigma,kappa,P);
    pChosen=P.choice1.*p1+(1-P.choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end

  function [p1,details]=probabilityAt(sigma,kappa,P)
    [p1,details]=maxpoolfit.choiceProbabilityDetectionGuessFromResponses( ...
      P.B,P.A,P.Rpre,cfg.adaptation,sigma,kappa,cfg,z1,z2, ...
      p.Results.GuessP_RDK1);
  end
end

function S=conditionDiagnostics(cache,mask,pRDK1,guessProbability)
T=cache.trials(mask,:); idx=find(mask);
vars=intersect({'domain','animal','probeOffsetDeg','stepMagnitudePC'}, ...
  T.Properties.VariableNames,'stable');
[groupID,keys]=findgroups(T(:,vars));
choice1=double(T.choiceRDK==1); change1=double(T.changeRDK==1);
p1=pRDK1(idx); pCorrect=change1.*p1+(1-change1).*(1-p1);
S=keys;
S.nTrials=splitapply(@numel,choice1,groupID);
S.observedAccuracy=splitapply(@mean,double(T.choiceRDK==T.changeRDK),groupID);
S.predictedAccuracy=splitapply(@mean,pCorrect,groupID);
S.accuracyResidual=S.observedAccuracy-S.predictedAccuracy;
S.observedPChooseRDK1=splitapply(@mean,choice1,groupID);
S.predictedPChooseRDK1=splitapply(@mean,p1,groupID);
S.meanGuessProbability=splitapply(@mean,guessProbability(idx),groupID);
end
