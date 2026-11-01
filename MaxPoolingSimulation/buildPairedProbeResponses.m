function model = buildPairedProbeResponses(bank,offsetDeg,yoked,cfg)
% BUILDPAIREDPROBERESPONSES  Construct response fields from coherence changes.
%
% Each mechanism receives its tuning-weighted mean step coherence minus its
% tuning-weighted mean prestep coherence.  BaseDriftCoherence,
% DriftPedestal, and ProbePedestal are therefore absent from the expressions
% below: they are present in both epochs and cancel exactly.  The Gaussian
% stream variables are scalar changes in epoch-mean coherence, not absolute
% step-epoch coherences or frame-resolved temporal noise.

streamDirections=[0,mod(offsetDeg,360),mod(-offsetDeg,360)];
delta=streamDirections(:)-bank.directionsDeg;
W=single(centeredMTSensitivity(delta,cfg.TuningSigmaDeg,cfg.CenterTuning));

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

changedDeltaStreams=single([cfg.DriftNoiseSD*bank.zDriftChanged, ...
    cfg.ProbeNoiseSD*zPlusChanged, ...
    cfg.ProbeNoiseSD*zMinusChanged]);
noChangeDeltaStreams=single([cfg.DriftNoiseSD*bank.zDriftNoChange, ...
    cfg.ProbeNoiseSD*zPlusNoChange, ...
    cfg.ProbeNoiseSD*zMinusNoChange]);

model=struct();
model.offsetDeg=offsetDeg;
model.yoked=yoked;
model.streamDirections=streamDirections;
model.W=W;
model.responseChangedBase=changedDeltaStreams*W+ ...
    single(cfg.InternalNoiseSD)*bank.internalChangedUnit;
model.responseNoChange=noChangeDeltaStreams*W+ ...
    single(cfg.InternalNoiseSD)*bank.internalNoChangeUnit;

driftDifference=cfg.DriftNoiseSD*(bank.zDriftChanged-bank.zDriftNoChange);
probe1Difference=cfg.ProbeNoiseSD*(zPlusChanged-zPlusNoChange);
if yoked
    model.predictors=[driftDifference,probe1Difference];
else
    probe2Difference=cfg.ProbeNoiseSD*(zMinusChanged-zMinusNoChange);
    model.predictors=[driftDifference,probe1Difference,probe2Difference];
end
end
