function model = buildDriftTemplateResponses(bank,offsetDeg,yoked, ...
    baseDriftCoherence,cancellationFraction,internalNoiseSD,cfg)
% BUILDDRIFTTEMPLATERESPONSES  Step landscapes after drift-template removal.

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

absoluteDrift=baseDriftCoherence+pcfg.DriftPedestal;
residualDrift=(1-cancellationFraction)*absoluteDrift;
residualMean=single([residualDrift,pcfg.ProbePedestal,pcfg.ProbePedestal])*W;
noiseScale=single(internalNoiseSD);

model=struct();
model.W=W;
model.streamDirections=streamDirections;
model.residualDriftCoherence=residualDrift;
model.cancelledDriftCoherence=cancellationFraction*absoluteDrift;
model.responseChangedBase=bsxfun(@plus,changedStreams*W,residualMean)+ ...
    noiseScale*bank.internalChangedUnit;
model.responseNoChange=bsxfun(@plus,noChangeStreams*W,residualMean)+ ...
    noiseScale*bank.internalNoChangeUnit;

driftDifference=pcfg.DriftNoiseSD*(bank.zDriftChanged-bank.zDriftNoChange);
probe1Difference=pcfg.ProbeNoiseSD*(plusChanged-plusNoChange);
if yoked
    model.predictors=[driftDifference,probe1Difference];
else
    probe2Difference=pcfg.ProbeNoiseSD*(minusChanged-minusNoChange);
    model.predictors=[driftDifference,probe1Difference,probe2Difference];
end
end
