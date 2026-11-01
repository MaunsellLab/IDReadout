function bank = makeSplitSignalTrialBank(nTrials,cfg,seed)
% MAKESPLITSIGNALTRIALBANK  Reusable noise for split-signal simulations.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
directionsDeg=0:cfg.BankSpacingDeg:(360-cfg.BankSpacingDeg);
nMechanisms=numel(directionsDeg);

bank=struct();
bank.nTrials=nTrials;
bank.directionsDeg=directionsDeg;
[bank.externalChangedUnit,bank.externalChangedFramesUnit]= ...
    integratedBinaryNoise(nTrials,2,cfg.NNoiseFrames);
[bank.externalNoChangeUnit,bank.externalNoChangeFramesUnit]= ...
    integratedBinaryNoise(nTrials,2,cfg.NNoiseFrames);
unitCfg=cfg; unitCfg.InternalNoiseSD=1;
bank.internalChangedUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.internalNoChangeUnit=generateCircularGPNoise(nTrials,nMechanisms,unitCfg);
bank.seed=seed;
end

function [z,frames]=integratedBinaryNoise(nTrials,nStreams,nFrames)
frames=zeros(nTrials,nStreams,nFrames,'single');
for frame=1:nFrames
    frames(:,:,frame)=single(2*(rand(nTrials,nStreams)>0.5)-1);
end
z=mean(frames,3);
end
