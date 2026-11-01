function joint = fitSplitSignalJointModels(bank,separationDeg,cfg)
% FITSPLITSIGNALJOINTMODELS  One slope/intercept per model across conditions.

y=[]; scores=struct('linear',[],'vector',[],'energy',[]);
softScores=cell(size(cfg.SoftmaxTemperatures));
for totalStep=cfg.TotalStepCoherences
    for allocation=cfg.AllocationFractions
        model=buildSplitSignalDesign(bank,separationDeg,allocation,totalStep, ...
            cfg.PrimaryStepNoiseHalfAmplitude, ...
            cfg.PrimaryPrestepPedestalFraction,cfg.PrimaryBalancePrestep, ...
            cfg.PrimaryAdaptationFraction,cfg);
        y=[y;model.choiceMax]; %#ok<AGROW>
        names=fieldnames(scores);
        for k=1:numel(names)
            name=names{k}; scores.(name)=[scores.(name);model.scores.(name)]; %#ok<AGROW>
        end
        for k=1:numel(cfg.SoftmaxTemperatures)
            temperature=cfg.SoftmaxTemperatures(k);
            score=splitSignalSoftMax(model.detChanged,temperature)- ...
                splitSignalSoftMax(model.detNoChange,temperature);
            softScores{k}=[softScores{k};score]; %#ok<AGROW>
        end
    end
end

names=fieldnames(scores); fits=struct();
for k=1:numel(names), name=names{k}; fits.(name)=fitLogisticIRLS(scores.(name),y); end
softFits=cell(size(softScores)); softLL=zeros(size(softScores));
for k=1:numel(softScores)
    softFits{k}=fitLogisticIRLS(softScores{k},y); softLL(k)=softFits{k}.logLikelihood;
end
[~,best]=max(softLL); fits.max=softFits{best};
joint=struct('separationDeg',separationDeg,'nTrials',numel(y),'fits',fits);
joint.bestSoftmaxTemperature=cfg.SoftmaxTemperatures(best);
joint.softmaxTemperatures=cfg.SoftmaxTemperatures;
joint.softmaxLogLikelihoods=softLL;
joint.deltaLL.maxMinusLinear=fits.max.logLikelihood-fits.linear.logLikelihood;
joint.deltaLL.maxMinusVector=fits.max.logLikelihood-fits.vector.logLikelihood;
joint.deltaLL.maxMinusEnergy=fits.max.logLikelihood-fits.energy.logLikelihood;
joint.pointDeltaLL.maxMinusLinear=fits.max.pointLogLikelihood-fits.linear.pointLogLikelihood;
joint.pointDeltaLL.maxMinusVector=fits.max.pointLogLikelihood-fits.vector.pointLogLikelihood;
joint.pointDeltaLL.maxMinusEnergy=fits.max.pointLogLikelihood-fits.energy.pointLogLikelihood;
end
