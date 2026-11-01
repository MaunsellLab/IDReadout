function bank = makeDensePairedTrialBank(nTrials,cfg,seed)
% MAKEDENSEPAIREDTRIALBANK  Common random inputs for all offsets/yoking modes.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
directionsDeg=0:cfg.BankSpacingDeg:(360-cfg.BankSpacingDeg);
nMechanisms=numel(directionsDeg);

bank=struct();
bank.nTrials=nTrials;
bank.directionsDeg=directionsDeg;
bank.zDriftChanged=randn(nTrials,1);
bank.zDriftNoChange=randn(nTrials,1);
bank.zProbe1Changed=randn(nTrials,1);
bank.zProbe1NoChange=randn(nTrials,1);
bank.zProbe2Changed=randn(nTrials,1);
bank.zProbe2NoChange=randn(nTrials,1);
unitNoiseCfg=cfg;
unitNoiseCfg.InternalNoiseSD=1;
bank.internalChangedUnit=generateCircularGPNoise(nTrials,nMechanisms,unitNoiseCfg);
bank.internalNoChangeUnit=generateCircularGPNoise(nTrials,nMechanisms,unitNoiseCfg);
bank.seed=seed;
end
