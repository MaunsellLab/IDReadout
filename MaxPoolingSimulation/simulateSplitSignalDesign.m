function result = simulateSplitSignalDesign(bank,separationDeg,allocation,totalStep, ...
    noiseHalfAmplitude,prestepFraction,balancePrestep,adaptationFraction,cfg)
% SIMULATESPLITSIGNALDESIGN  Performance, gains, and alternative fits.

model=buildSplitSignalDesign(bank,separationDeg,allocation,totalStep, ...
    noiseHalfAmplitude,prestepFraction,balancePrestep,adaptationFraction,cfg);
names={'max','linear','vector','energy'};
fits=struct();
for k=1:numel(names)
    name=names{k}; fits.(name)=fitLogisticIRLS(model.scores.(name),model.choiceMax);
end
gainFit=fitProbitIRLS(model.noiseDifference,model.choiceMax);
if allocation>=0.5
    winnerIndex=1; loserIndex=2;
else
    winnerIndex=2; loserIndex=1;
end
winnerGain=gainFit.beta(winnerIndex); loserGain=gainFit.beta(loserIndex);

result=struct();
result.separationDeg=separationDeg;
result.allocation=allocation;
result.totalStep=totalStep;
result.noiseHalfAmplitude=noiseHalfAmplitude;
result.prestepFraction=prestepFraction;
result.balancePrestep=balancePrestep;
result.adaptationFraction=adaptationFraction;
result.performance=mean(model.choiceMax);
result.meanStreams=model.meanStreams;
result.prestepPedestal=model.prestepPedestal;
result.commonPedestalRise=model.commonPedestalRise;
result.allocatableSignal=model.allocatableSignal;
result.balanceDirection=model.balanceDirection;
result.balanceCoherence=model.balanceCoherence;
result.gainFit=gainFit;
result.winnerGain=winnerGain;
result.loserGain=loserGain;
result.loserWinnerGainRatio=loserGain/winnerGain;
compactFits=fits;
for k=1:numel(names)
    name=names{k};
    compactFits.(name)=rmfield(compactFits.(name),{'probability','pointLogLikelihood'});
end
result.modelFits=compactFits;
result.deltaLL.maxMinusLinear=fits.max.logLikelihood-fits.linear.logLikelihood;
result.deltaLL.maxMinusVector=fits.max.logLikelihood-fits.vector.logLikelihood;
result.deltaLL.maxMinusEnergy=fits.max.logLikelihood-fits.energy.logLikelihood;
end
