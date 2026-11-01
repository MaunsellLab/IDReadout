function model = preparePoissonPairedProbeCondition(bank,offsetDeg,yoked, ...
    stepCoherence,baseDriftCoherence,cfg)
% PREPAREPOISSONPAIREDPROBECONDITION  Mean delta and N=1 Poisson noise.
%
% Expected absolute rates are constructed for both epochs.  Conditional count
% noise is approximated by correlated Gaussian fields with the exact pooled-
% Poisson marginal variance.  model.noise*Unit fields correspond to one
% independent-neuron equivalent and are divided by sqrt(EffectivePoolSize)
% when choices are simulated.

pcfg=cfg.PairedConfig;
streamDirections=[0,mod(offsetDeg,360),mod(-offsetDeg,360)];
delta=streamDirections(:)-bank.directionsDeg;
W=single(centeredMTSensitivity(delta,pcfg.TuningSigmaDeg,pcfg.CenterTuning));

if yoked
    zPlusChanged=bank.zProbe1Changed;
    zMinusChanged=zPlusChanged;
    zPlusNoChange=bank.zProbe1NoChange;
    zMinusNoChange=zPlusNoChange;
else
    zPlusChanged=bank.zProbe1Changed;
    zMinusChanged=bank.zProbe2Changed;
    zPlusNoChange=bank.zProbe1NoChange;
    zMinusNoChange=bank.zProbe2NoChange;
end

changedDeltaStreams=single([pcfg.DriftNoiseSD*bank.zDriftChanged, ...
    pcfg.ProbeNoiseSD*zPlusChanged,pcfg.ProbeNoiseSD*zMinusChanged]);
noChangeDeltaStreams=single([pcfg.DriftNoiseSD*bank.zDriftNoChange, ...
    pcfg.ProbeNoiseSD*zPlusNoChange,pcfg.ProbeNoiseSD*zMinusNoChange]);

alignedSensitivity=double(W(1,1));
rateScale=single(cfg.AlignedRateSlopeSpPerSecPerPct/alignedSensitivity);
preStreams=single([baseDriftCoherence+pcfg.DriftPedestal, ...
    pcfg.ProbePedestal,pcfg.ProbePedestal]);
preModel=preStreams*W;
preRateRaw=single(cfg.RandomMotionRateSpPerSec)+rateScale*preModel;

changedMeanModel=changedDeltaStreams*W+single(stepCoherence)*W(1,:);
noChangeMeanModel=noChangeDeltaStreams*W;
changedStepRateRaw=bsxfun(@plus,preRateRaw,rateScale*changedMeanModel);
noChangeStepRateRaw=bsxfun(@plus,preRateRaw,rateScale*noChangeMeanModel);

minimumRate=single(cfg.MinimumRateSpPerSec);
preRate=max(preRateRaw,minimumRate);
changedStepRate=max(changedStepRateRaw,minimumRate);
noChangeStepRate=max(noChangeStepRateRaw,minimumRate);

% Use the rates after the positivity bound for both the mean and variance.
meanChanged=bsxfun(@minus,changedStepRate,preRate)/rateScale;
meanNoChange=bsxfun(@minus,noChangeStepRate,preRate)/rateScale;

preSD=sqrt(single(cfg.FanoFactor)*preRate/single(cfg.PreStepDurationSec))/rateScale;
changedStepSD=sqrt(single(cfg.FanoFactor)*changedStepRate/ ...
    single(cfg.StepDurationSec))/rateScale;
noChangeStepSD=sqrt(single(cfg.FanoFactor)*noChangeStepRate/ ...
    single(cfg.StepDurationSec))/rateScale;

noiseChangedUnit=changedStepSD.*bank.stepChangedUnit- ...
    bsxfun(@times,bank.preChangedUnit,preSD);
noiseNoChangeUnit=noChangeStepSD.*bank.stepNoChangeUnit- ...
    bsxfun(@times,bank.preNoChangeUnit,preSD);

driftDifference=pcfg.DriftNoiseSD*(bank.zDriftChanged-bank.zDriftNoChange);
probe1Difference=pcfg.ProbeNoiseSD*(zPlusChanged-zPlusNoChange);
if yoked
    predictors=[driftDifference,probe1Difference];
else
    probe2Difference=pcfg.ProbeNoiseSD*(zMinusChanged-zMinusNoChange);
    predictors=[driftDifference,probe1Difference,probe2Difference];
end

[~,driftIndex]=min(abs(bank.directionsDeg));
[~,probeIndex]=min(abs(mod(bank.directionsDeg-offsetDeg+180,360)-180));
model=struct();
model.offsetDeg=offsetDeg;
model.yoked=yoked;
model.streamDirections=streamDirections;
model.W=W;
model.meanChanged=meanChanged;
model.meanNoChange=meanNoChange;
model.noiseChangedUnit=noiseChangedUnit;
model.noiseNoChangeUnit=noiseNoChangeUnit;
model.predictors=predictors;
model.preRateDrift=double(preRate(driftIndex));
model.preRateProbe=double(preRate(probeIndex));
model.preRateMinimum=double(min(preRate));
model.preRateMaximum=double(max(preRate));
model.changedClampFraction=mean(changedStepRateRaw(:)<minimumRate);
model.noChangeClampFraction=mean(noChangeStepRateRaw(:)<minimumRate);
model.preClampFraction=mean(preRateRaw(:)<minimumRate);
model.rateScaleSpPerSecPerModelUnit=double(rateScale);
end
