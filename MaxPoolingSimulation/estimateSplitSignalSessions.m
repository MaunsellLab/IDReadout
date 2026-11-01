function power = estimateSplitSignalSessions(pointDeltaLL,cfg,seed)
% ESTIMATESPLITSIGNALSESSIONS  Bootstrap session counts from pointwise delta LL.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
pointDeltaLL=double(pointDeltaLL(:)); n=numel(pointDeltaLL);
sessions=(1:cfg.MaxSessions)'; probability=zeros(size(sessions));
running=zeros(cfg.NPowerReplicates,1);
for s=1:numel(sessions)
    for rep=1:cfg.NPowerReplicates
        indices=randi(n,cfg.TrialsPerSession,1);
        running(rep)=running(rep)+sum(pointDeltaLL(indices));
    end
    probability(s)=mean(running>=cfg.PowerCriterionDeltaLL);
end
first=find(probability>=cfg.TargetPower,1,'first');
if isempty(first), required=NaN; else, required=sessions(first); end
power=struct('sessions',sessions,'probability',probability, ...
    'requiredSessions',required,'criterionDeltaLL',cfg.PowerCriterionDeltaLL, ...
    'targetPower',cfg.TargetPower);
end
