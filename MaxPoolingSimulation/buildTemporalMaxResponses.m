function model = buildTemporalMaxResponses(bank,offsetDeg,yoked, ...
    stepCoherence,baseDriftCoherence,cfg)
% BUILDTEMPORALMAXRESPONSES  Absolute epoch landscapes before separate maxima.

pcfg=cfg.PairedConfig;
streamDirections=[0,mod(offsetDeg,360),mod(-offsetDeg,360)];
delta=streamDirections(:)-bank.directionsDeg;
W=single(centeredMTSensitivity(delta,pcfg.TuningSigmaDeg,pcfg.CenterTuning));

if yoked
    plusChanged=bank.zProbe1Changed; minusChanged=plusChanged;
    plusNoChange=bank.zProbe1NoChange; minusNoChange=plusNoChange;
else
    plusChanged=bank.zProbe1Changed; minusChanged=bank.zProbe2Changed;
    plusNoChange=bank.zProbe1NoChange; minusNoChange=bank.zProbe2NoChange;
end
changedStreams=single([pcfg.DriftNoiseSD*bank.zDriftChanged, ...
    pcfg.ProbeNoiseSD*plusChanged,pcfg.ProbeNoiseSD*minusChanged]);
noChangeStreams=single([pcfg.DriftNoiseSD*bank.zDriftNoChange, ...
    pcfg.ProbeNoiseSD*plusNoChange,pcfg.ProbeNoiseSD*minusNoChange]);

preStreams=single([baseDriftCoherence+pcfg.DriftPedestal, ...
    pcfg.ProbePedestal,pcfg.ProbePedestal]);
preMean=preStreams*W;
model=struct();
model.W=W;
model.streamDirections=streamDirections;
model.preMean=preMean;
model.stepChangedMean=bsxfun(@plus,changedStreams*W,preMean)+ ...
    single(stepCoherence)*W(1,:);
model.stepNoChangeMean=bsxfun(@plus,noChangeStreams*W,preMean);
preScale=single(sqrt(cfg.StepDurationSec/cfg.PreStepDurationSec));
model.preChangedUnit=preScale*bank.preChangedUnit;
model.preNoChangeUnit=preScale*bank.preNoChangeUnit;
model.stepChangedUnit=bank.stepChangedUnit;
model.stepNoChangeUnit=bank.stepNoChangeUnit;

driftDifference=pcfg.DriftNoiseSD*(bank.zDriftChanged-bank.zDriftNoChange);
probe1Difference=pcfg.ProbeNoiseSD*(plusChanged-plusNoChange);
if yoked
    model.predictors=[driftDifference,probe1Difference];
else
    probe2Difference=pcfg.ProbeNoiseSD*(minusChanged-minusNoChange);
    model.predictors=[driftDifference,probe1Difference,probe2Difference];
end
end
