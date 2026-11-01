function gains = estimateSplitSignalRecoveryGains(design,trueProbability,cfg,seed)
% ESTIMATESPLITSIGNALRECOVERYGAINS  Frame-resolved gains with empirical SEs.

oldState=rng; cleanup=onCleanup(@()rng(oldState)); rng(seed,'twister');
nAlloc=numel(cfg.AllocationFractions); nFrames=cfg.NNoiseFrames;
g1=nan(cfg.NGainReplicates,nAlloc,nFrames);
g2=nan(cfg.NGainReplicates,nAlloc,nFrames);
for rep=1:cfg.NGainReplicates
    outcome=rand(design.nTrials,1)<trueProbability;
    for a=1:nAlloc
        index=find(abs(design.allocation-cfg.AllocationFractions(a))<1e-9);
        f1=reshape(design.frameDifference(index,1,:),numel(index),nFrames);
        f2=reshape(design.frameDifference(index,2,:),numel(index),nFrames);
        predictors=[design.totalStep(index),f1,f2];
        fit=fitProbitIRLS(predictors,outcome(index));
        g1(rep,a,:)=reshape(fit.beta(2:(1+nFrames)),1,1,nFrames);
        g2(rep,a,:)=reshape(fit.beta((2+nFrames):(1+2*nFrames)),1,1,nFrames);
    end
end
gains=struct(); gains.allocationFractions=cfg.AllocationFractions;
gains.stream1FrameMean=squeeze(mean(g1,1));
gains.stream2FrameMean=squeeze(mean(g2,1));
gains.stream1FrameSE=squeeze(std(g1,0,1));
gains.stream2FrameSE=squeeze(std(g2,0,1));
gains.stream1Mean=mean(gains.stream1FrameMean,2);
gains.stream2Mean=mean(gains.stream2FrameMean,2);
gains.stream1MeanSE=std(mean(g1,3),0,1)';
gains.stream2MeanSE=std(mean(g2,3),0,1)';
gains.equivalentSessions=design.nTrials/cfg.TrialsPerSession;
end
