function result=simulateUpstreamChoiceCondition(condition,bank,noiseModel, ...
  mtModel,noiseSDByTrial,nReplicates,stream,varargin)
% simulateUpstreamChoiceCondition  Replicated dense-bank hard-max choices.

p=inputParser;
addParameter(p,'RectifyCandidates',true,@(x) islogical(x)&&isscalar(x));
addParameter(p,'ChunkSize',3000,@(x) isnumeric(x)&&isscalar(x)&&x>=1&&fix(x)==x);
parse(p,varargin{:});opts=p.Results;

nTrials=numel(condition.trialIndex);
noiseSDByTrial=double(noiseSDByTrial(:));
if numel(noiseSDByTrial)~=nTrials||any(~isfinite(noiseSDByTrial))|| ...
    any(noiseSDByTrial<0)
  error('simulateUpstreamChoiceCondition:BadNoiseScale', ...
    'noiseSDByTrial must contain one finite nonnegative value per trial.');
end
if ~(isscalar(nReplicates)&&nReplicates>=1&&fix(nReplicates)==nReplicates)
  error('simulateUpstreamChoiceCondition:BadReplicates', ...
    'nReplicates must be a positive integer.');
end

C=evaluateDenseMTStepCandidates(condition.componentDirectionsDeg, ...
  condition.changeComponentCoherence,bank,mtModel,'RectifyCandidates',false);
N=evaluateDenseMTStepCandidates(condition.componentDirectionsDeg, ...
  condition.noChangeComponentCoherence,bank,mtModel,'RectifyCandidates',false);
nCandidates=numel(bank.candidateDirectionsDeg);
correctReplicates=false(nTrials,nReplicates);
replicateCorrectFraction=nan(nReplicates,1);
changeWinnerCounts=zeros(nCandidates,1);
noChangeWinnerCounts=zeros(nCandidates,1);

for iRep=1:nReplicates
  for first=1:opts.ChunkSize:nTrials
    rows=first:min(first+opts.ChunkSize-1,nTrials);
    n=numel(rows);scale=noiseSDByTrial(rows);
    z=randn(stream,n,noiseModel.rank)*noiseModel.factor';
    changeActivation=C.evidence(rows,:)+scale.*z;
    z=randn(stream,n,noiseModel.rank)*noiseModel.factor';
    noChangeActivation=N.evidence(rows,:)+scale.*z;
    if opts.RectifyCandidates
      changeActivation=max(changeActivation,0);
      noChangeActivation=max(noChangeActivation,0);
    end
    [changeMax,changeWinner]=max(changeActivation,[],2);
    [noChangeMax,noChangeWinner]=max(noChangeActivation,[],2);
    correctReplicates(rows,iRep)=changeMax>noChangeMax;
    changeWinnerCounts=changeWinnerCounts+accumarray( ...
      changeWinner,1,[nCandidates 1]);
    noChangeWinnerCounts=noChangeWinnerCounts+accumarray( ...
      noChangeWinner,1,[nCandidates 1]);
  end
  replicateCorrectFraction(iRep)=mean(correctReplicates(:,iRep));
end

result=struct();
result.dataset=condition.dataset;result.offsetDeg=condition.offsetDeg;
result.sessionIndex=condition.sessionIndex;result.trialIndex=condition.trialIndex;
result.experimentalCorrect=condition.experimentalCorrect;
result.correctReplicates=correctReplicates;
result.trialProbabilityCorrect=mean(correctReplicates,2);
result.replicateCorrectFraction=replicateCorrectFraction;
result.candidateDirectionsDeg=bank.candidateDirectionsDeg;
result.changeWinnerFraction=changeWinnerCounts./(nTrials*nReplicates);
result.noChangeWinnerFraction=noChangeWinnerCounts./(nTrials*nReplicates);
result.rectifyCandidates=opts.RectifyCandidates;
result.poolingMode="hardMax";
end
