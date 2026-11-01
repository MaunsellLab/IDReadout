function design = makeSplitSignalRecoveryDesign(thresholds,cfg,seed)
% MAKESPLITSIGNALRECOVERYDESIGN  Balanced trial bank around model thresholds.

nConditions=numel(cfg.AllocationFractions)*numel(cfg.SignalLevelMultipliers);
nPer=ceil(cfg.NStimulusTrials/nConditions); nTotal=nPer*nConditions;
bank=makeSplitSignalTrialBank(nTotal,cfg,seed);
nMechanisms=numel(bank.directionsDeg);
detChanged=zeros(nTotal,nMechanisms,'single');
detNoChange=zeros(nTotal,nMechanisms,'single');
frameDifference=zeros(nTotal,2,cfg.NNoiseFrames,'single');
allocation=zeros(nTotal,1); totalStep=zeros(nTotal,1); level=zeros(nTotal,1);

first=1;
for a=1:numel(cfg.AllocationFractions)
    for l=1:numel(cfg.SignalLevelMultipliers)
        index=first:(first+nPer-1); first=first+nPer;
        subBank=subsetSplitSignalTrialBank(bank,index);
        step=thresholds(a)*cfg.SignalLevelMultipliers(l);
        model=buildSplitSignalDesign(subBank,cfg.SeparationDeg, ...
            cfg.AllocationFractions(a),step,cfg.StepNoiseHalfAmplitude, ...
            cfg.PrestepPedestalFraction,cfg.BalancePrestep, ...
            cfg.AdaptationFraction,cfg);
        detChanged(index,:)=model.detChanged;
        detNoChange(index,:)=model.detNoChange;
        frameDifference(index,:,:)=single(cfg.StepNoiseHalfAmplitude)* ...
            (subBank.externalChangedFramesUnit-subBank.externalNoChangeFramesUnit);
        allocation(index)=cfg.AllocationFractions(a);
        totalStep(index)=step; level(index)=cfg.SignalLevelMultipliers(l);
    end
end

design=struct('nTrials',nTotal,'directionsDeg',bank.directionsDeg, ...
    'detChanged',detChanged,'detNoChange',detNoChange, ...
    'frameDifference',frameDifference,'allocation',allocation, ...
    'totalStep',totalStep,'levelMultiplier',level,'thresholds',thresholds);
end
