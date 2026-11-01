function recovered = recoverMTReadoutPNormGains(result,summaries)
%RECOVERMTREADOUTPNORMGAINS Use the EXISTING IDS predictors and gain fitter.
% summaries: output from plotIDSSummary('Animal','Meetz','IncludeIDQ',false)
% or acrossSummary structs loaded from its saved output files. The frozen
% trialTable contains the exact existing looStep (or other selected)
% predictors and the observed psychometric calibration.
% Requires the user's fitIDSNoiseGain, idqWeibullP, and
% idqWeibullAlphaForThreshold on the MATLAB path (plus their toolboxes).
% Replacing Bernoulli outcomes by model p(correct) gives the expected
% Bernoulli log likelihood; it avoids extra simulated-choice sampling noise.
assert(exist('fitIDSNoiseGain','file')==2,'Add the current IDS gain-analysis code to the path.');
if isstruct(summaries), summaries=num2cell(summaries); end
modelT=result.trials;
modelKeys=string(modelT.sessionID)+"/"+string(modelT.parentTrialIdx);
assert(numel(unique(modelKeys))==height(modelT),'Duplicate model trial keys.');
recovered=cell(numel(summaries),1);
for s=1:numel(summaries)
  S=summaries{s}; T=S.trialTable;
  assert(strcmpi(string(S.animal),result.config.animal),'Animal mismatch.');
  if isfield(S,'sessionRecords') && ismember('sourceDomain',S.sessionRecords.Properties.VariableNames)
    assert(all(strcmpi(string(S.sessionRecords.sourceDomain),'IDS')), ...
      'Recreate summaries with IncludeIDQ=false for this initial IDS model.');
  end
  required={'sourceSessionID','parentTrialIdx','noisePredDrift','noisePredPlus','noisePredMinus'};
  assert(all(ismember(required,T.Properties.VariableNames)),'IDS summary lacks frozen predictor fields.');
  keys=string(T.sourceSessionID)+"/"+string(T.parentTrialIdx);
  [found,idx]=ismember(keys,modelKeys);
  assert(all(found),'%d summary trials lack model predictions. Supply matching full sessions.',nnz(~found));
  assert(all(abs(double(T.stepCoh)-modelT.stepMagnitudePC(idx))<1e-4), ...
    'Step magnitudes differ: check corrected IDS descriptors.');
  assert(all(double(T.sideIndex)==modelT.changeRDK(idx)), 'Changed-side coding differs.');
  assert(all(logical(T.correct)==modelT.correct(idx)), 'Observed outcomes differ.');
  assert(all(modelT.isIncrement(idx)), ...
    'Initial gain recovery is for increment trials; supply increment-only summaries.');
  fits=cell(height(result.fits),1);
  for k=1:height(result.fits)
    M=T; M.correct=result.predictions{k}.pCorrect(idx);
    fits{k}=fitIDSNoiseGain(M,S.psychFit.sessionFits,S.psychFit.alignedWeibull, ...
      S.psychFit.targetPerformance);
  end
  observed=fitIDSNoiseGain(T,S.psychFit.sessionFits,S.psychFit.alignedWeibull, ...
    S.psychFit.targetPerformance);
  recovered{s}=struct('probeOffsetDeg',S.probeOffsetDeg,'observed',observed, ...
    'predicted',{fits},'p',result.fits.p,'noisePredictorType',S.noisePredictorType, ...
    'nTrials',height(T),'nSessions',numel(unique(string(T.sourceSessionID))), ...
    'psychometricCalibration',"fixed observed calibration", ...
    'predictorContract',"existing frozen IDS predictors, unchanged");
end
end
