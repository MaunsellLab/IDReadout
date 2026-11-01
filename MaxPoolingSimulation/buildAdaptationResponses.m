function model = buildAdaptationResponses(bank,offsetDeg,yoked, ...
    baseDriftCoherence,aftereffectAt60,internalNoiseSD,cfg)
% BUILDADAPTATIONRESPONSES  Step landscapes after cancellation and adaptation.

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
actualAftereffect=aftereffectAt60*absoluteDrift/cfg.AdaptationReferenceCoherence;
adaptationDistance=mod(double(bank.directionsDeg)+180,360)-180;
adaptationProfile=single(exp(-0.5*(adaptationDistance/ ...
    cfg.AdaptationWidthDeg).^2));

% Exact perceptual cancellation removes the sustained drift signal.  Probe
% pedestals remain.  Adaptation is a nonnegative suppressive offset: quiet
% far-direction mechanisms are spared rather than literally excited.
residualMean=single([0,pcfg.ProbePedestal,pcfg.ProbePedestal])*W- ...
    single(actualAftereffect)*adaptationProfile;
noiseScale=single(internalNoiseSD);

model=struct();
model.W=W;
model.streamDirections=streamDirections;
model.adaptationProfile=adaptationProfile;
model.actualAftereffect=actualAftereffect;
model.cancelledDriftCoherence=absoluteDrift;
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
