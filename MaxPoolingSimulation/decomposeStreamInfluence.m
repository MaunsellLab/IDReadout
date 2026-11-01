function d = decomposeStreamInfluence(responseChanged,responseNoChange,profile,delta)
% DECOMPOSESTREAMINFLUENCE  Fixed-winner and winner-switch evidence effects.

[oldChanged,winnerChanged]=max(responseChanged,[],2);
[oldNoChange,winnerNoChange]=max(responseNoChange,[],2);
fixedChanged=mean(double(profile(winnerChanged)));
fixedNoChange=mean(double(profile(winnerNoChange)));

[newChanged,newWinnerChanged]=max(bsxfun(@plus,responseChanged, ...
    single(delta)*profile),[],2);
[newNoChange,newWinnerNoChange]=max(bsxfun(@plus,responseNoChange, ...
    single(delta)*profile),[],2);
totalChanged=mean(double(newChanged-oldChanged))/delta;
totalNoChange=mean(double(newNoChange-oldNoChange))/delta;

d=struct();
d.changedFixed=fixedChanged;
d.changedSwitch=totalChanged-fixedChanged;
d.noChangeFixed=fixedNoChange;
d.noChangeSwitch=totalNoChange-fixedNoChange;
d.decisionFixed=0.5*(fixedChanged+fixedNoChange);
d.decisionSwitch=0.5*(d.changedSwitch+d.noChangeSwitch);
d.decisionTotal=d.decisionFixed+d.decisionSwitch;
d.changedSwitchProbability=mean(newWinnerChanged~=winnerChanged);
d.noChangeSwitchProbability=mean(newWinnerNoChange~=winnerNoChange);
end

