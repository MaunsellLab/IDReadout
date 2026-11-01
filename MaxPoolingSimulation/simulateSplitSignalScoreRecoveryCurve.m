function recovery = simulateSplitSignalScoreRecoveryCurve(trueProbability, ...
    scores,trueModelIndex,cfg,seed)
% SIMULATESPLITSIGNALSCORESRECOVERYCURVE  Refit all score models each time.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
n=numel(trueProbability); nModels=numel(cfg.PoolingModels);
winnerCounts=zeros(cfg.MaxRecoverySessions,nModels);
for rep=1:cfg.NRecoveryReplicates
    allIndex=[]; allOutcome=[];
    for session=1:cfg.MaxRecoverySessions
        index=randi(n,cfg.TrialsPerSession,1);
        outcome=rand(cfg.TrialsPerSession,1)<trueProbability(index);
        allIndex=[allIndex;index]; allOutcome=[allOutcome;outcome]; %#ok<AGROW>
        subScores=struct('max',scores.max(allIndex,:), ...
            'vector',scores.vector(allIndex),'energy',scores.energy(allIndex));
        fit=fitSplitSignalScoreModels(allOutcome,subScores,cfg);
        winnerCounts(session,fit.winnerIndex)=winnerCounts(session,fit.winnerIndex)+1;
    end
end
probability=winnerCounts/cfg.NRecoveryReplicates;
diagonal=probability(:,trueModelIndex);
first=find(diagonal>=cfg.TargetRecoveryProbability,1,'first');
if isempty(first), required=NaN; else, required=first; end
recovery=struct('trueModelIndex',trueModelIndex,'probability',probability, ...
    'diagonalProbability',diagonal,'requiredSessions',required);
end
