function d = decomposeTemporalMaxInfluence(preWinnerChanged,preWinnerNoChange, ...
    stepResponseChanged,stepResponseNoChange,profile,delta)
% DECOMPOSETEMPORALMAXINFLUENCE  Same/different temporal winners and switches.

[oldC,winnerC]=max(stepResponseChanged,[],2);
[oldN,winnerN]=max(stepResponseNoChange,[],2);
fixedC=double(profile(winnerC)); fixedC=fixedC(:);
fixedN=double(profile(winnerN)); fixedN=fixedN(:);
sameC=winnerC==preWinnerChanged;
sameN=winnerN==preWinnerNoChange;
[newC,newWinnerC]=max(bsxfun(@plus,stepResponseChanged, ...
    single(delta)*profile),[],2);
[newN,newWinnerN]=max(bsxfun(@plus,stepResponseNoChange, ...
    single(delta)*profile),[],2);
totalC=mean(double(newC-oldC))/delta;
totalN=mean(double(newN-oldN))/delta;

d=struct();
d.changedTemporalSameFixed=mean(fixedC.*sameC);
d.changedTemporalDifferentFixed=mean(fixedC.*~sameC);
d.changedPerturbationSwitch=totalC-mean(fixedC);
d.noChangeTemporalSameFixed=mean(fixedN.*sameN);
d.noChangeTemporalDifferentFixed=mean(fixedN.*~sameN);
d.noChangePerturbationSwitch=totalN-mean(fixedN);
d.decisionTemporalSameFixed=0.5*(d.changedTemporalSameFixed+ ...
    d.noChangeTemporalSameFixed);
d.decisionTemporalDifferentFixed=0.5*(d.changedTemporalDifferentFixed+ ...
    d.noChangeTemporalDifferentFixed);
d.decisionPerturbationSwitch=0.5*(d.changedPerturbationSwitch+ ...
    d.noChangePerturbationSwitch);
d.decisionTotal=d.decisionTemporalSameFixed+ ...
    d.decisionTemporalDifferentFixed+d.decisionPerturbationSwitch;
d.changedTemporalSameProbability=mean(sameC);
d.noChangeTemporalSameProbability=mean(sameN);
d.changedPerturbationSwitchProbability=mean(newWinnerC~=winnerC);
d.noChangePerturbationSwitchProbability=mean(newWinnerN~=winnerN);
end
