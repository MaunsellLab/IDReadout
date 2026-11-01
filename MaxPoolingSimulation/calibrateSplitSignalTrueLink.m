function link = calibrateSplitSignalTrueLink(thresholds,trueModelIndex,cfg,seed)
% CALIBRATESPLITSIGNALTRUELINK  Smooth link learned on independent simulations.

linkCfg=cfg; linkCfg.NStimulusTrials=cfg.NLinkCalibrationTrials;
design=makeSplitSignalRecoveryDesign(thresholds,linkCfg,seed);
unitCfg=cfg; unitCfg.InternalNoiseSD=1;
unitChanged=generateCircularGPNoise(design.nTrials,numel(design.directionsDeg),unitCfg);
unitNoChange=generateCircularGPNoise(design.nTrials,numel(design.directionsDeg),unitCfg);
changed=design.detChanged+single(cfg.TrueInternalNoiseSD)*unitChanged;
noChange=design.detNoChange+single(cfg.TrueInternalNoiseSD)*unitNoChange;
name=cfg.PoolingModels{trueModelIndex};
choice=poolSplitSignalResponse(changed,design.directionsDeg,name)> ...
    poolSplitSignalResponse(noChange,design.directionsDeg,name);
scores=computeSplitSignalScoreSet(design,cfg);
allFits=fitSplitSignalScoreModels(choice,scores,cfg);
link=struct('modelIndex',trueModelIndex,'modelName',name, ...
    'fit',allFits.detail{trueModelIndex},'bestTemperature',NaN);
if strcmpi(name,'max'), link.bestTemperature=allFits.bestTemperature; end
end
