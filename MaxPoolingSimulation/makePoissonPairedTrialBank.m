function bank = makePoissonPairedTrialBank(nTrials,cfg,seed)
% MAKEPOISSONPAIREDTRIALBANK  Common stream draws and four smooth count fields.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
pcfg=cfg.PairedConfig;
directionsDeg=0:pcfg.BankSpacingDeg:(360-pcfg.BankSpacingDeg);
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

unitCfg=pcfg;
unitCfg.InternalNoiseSD=1;
bank.preChangedUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.stepChangedUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.preNoChangeUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.stepNoChangeUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.seed=seed;
end
