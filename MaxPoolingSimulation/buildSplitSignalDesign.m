function model = buildSplitSignalDesign(bank,separationDeg,allocation,totalStep, ...
    noiseHalfAmplitude,prestepFraction,balancePrestep,adaptationFraction,cfg)
% BUILDSPLITSIGNALDESIGN  Responses and model scores for one condition.
%
% Carrier directions are 0 and separationDeg.  Both receive a positive mean
% step on every signal trial.  The prestep pedestal is a fraction of the
% step-noise half-amplitude.  Optional balancing motion opposes the vector
% sum of the two prestep pedestals and remains constant during the step.

carrierDirections=[0 mod(separationDeg,360)];
delta=carrierDirections(:)-bank.directionsDeg;
H=single(centeredMTSensitivity(delta,cfg.TuningSigmaDeg,cfg.CenterTuning));
T=single(centeredMTSensitivity(delta,cfg.TuningSigmaDeg,false));

prestepPedestal=prestepFraction*noiseHalfAmplitude;
preResponse=single(prestepPedestal)*sum(T,1);
balanceDirection=NaN; balanceCoherence=0;
if balancePrestep && prestepPedestal>0
    resultant=prestepPedestal*(exp(1i*deg2rad(carrierDirections(1)))+ ...
        exp(1i*deg2rad(carrierDirections(2))));
    balanceCoherence=abs(resultant);
    balanceDirection=mod(rad2deg(angle(resultant))+180,360);
    balanceDelta=balanceDirection-bank.directionsDeg;
    balanceT=single(centeredMTSensitivity(balanceDelta, ...
        cfg.TuningSigmaDeg,false));
    preResponse=preResponse+single(balanceCoherence)*balanceT;
end

adaptation=single(adaptationFraction)*preResponse;
residualPre=-single(cfg.TemplateCancellationFraction)*adaptation;
commonPedestalRise=(1-prestepFraction)*noiseHalfAmplitude;
allocatableSignal=totalStep-2*commonPedestalRise;
if allocatableSignal<0
    error('Total step is smaller than the two required pedestal rises.');
end
meanStreams=single(commonPedestalRise+ ...
    allocatableSignal*[allocation 1-allocation]);
commonNoChange=single([commonPedestalRise commonPedestalRise]);
changedStreams=bsxfun(@plus, ...
    single(noiseHalfAmplitude)*bank.externalChangedUnit,meanStreams);
noChangeStreams=bsxfun(@plus, ...
    single(noiseHalfAmplitude)*bank.externalNoChangeUnit,commonNoChange);

detChanged=bsxfun(@plus,changedStreams*H,residualPre);
detNoChange=bsxfun(@plus,noChangeStreams*H,residualPre);
responseChanged=detChanged+single(cfg.InternalNoiseSD)*bank.internalChangedUnit;
responseNoChange=detNoChange+single(cfg.InternalNoiseSD)*bank.internalNoChangeUnit;

maxChanged=max(responseChanged,[],2);
maxNoChange=max(responseNoChange,[],2);
choiceMax=maxChanged>maxNoChange;

zChanged=changedStreams*single(exp(1i*deg2rad(carrierDirections(:))));
zNoChange=noChangeStreams*single(exp(1i*deg2rad(carrierDirections(:))));
scoreVector=abs(zChanged)-abs(zNoChange);
scoreEnergy=sum(changedStreams.^2,2)-sum(noChangeStreams.^2,2);
scoreLinear=sum(changedStreams,2)-sum(noChangeStreams,2);
scoreMax=splitSignalSoftMax(detChanged,cfg.PrimarySoftmaxTemperature)- ...
    splitSignalSoftMax(detNoChange,cfg.PrimarySoftmaxTemperature);

model=struct();
model.choiceMax=choiceMax;
model.correct=choiceMax;
model.maxDecision=maxChanged-maxNoChange;
model.scores=struct('max',scoreMax,'linear',scoreLinear, ...
    'vector',scoreVector,'energy',scoreEnergy);
model.noiseDifference=single(noiseHalfAmplitude)* ...
    (bank.externalChangedUnit-bank.externalNoChangeUnit);
model.meanStreams=meanStreams;
model.prestepPedestal=prestepPedestal;
model.commonPedestalRise=commonPedestalRise;
model.allocatableSignal=allocatableSignal;
model.balanceDirection=balanceDirection;
model.balanceCoherence=balanceCoherence;
model.adaptation=adaptation;
model.H=H;
model.detChanged=detChanged;
model.detNoChange=detNoChange;
end
