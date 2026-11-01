function diagnostic = decisionRuleGainDiagnostic(cache,primaryResult,varargin)
%DECISIONRULEGAINDIAGNOSTIC Compare absolute-step and directional-change gains.

p=inputParser;
p.addParameter('Domain',"IDR",@(x)isstring(x)||ischar(x));
p.addParameter('Animal',"Meetz",@(x)isstring(x)||ischar(x));
p.addParameter('ProbeOffsetDeg',25,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('StepMagnitudePC',23,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('FitDraws',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=2));
p.addParameter('PredictionDraws',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)&&x>=2));
p.addParameter('MinimumTrials',500,@(x)isnumeric(x)&&isscalar(x)&&x>=20);
p.addParameter('MinimumSessions',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(isfield(primaryResult,'finalFit')&&isfield(primaryResult,'prediction'), ...
  'Primary result must contain finalFit and prediction.');
primaryFit=primaryResult.finalFit;
fitDraws=p.Results.FitDraws;
if isempty(fitDraws), fitDraws=primaryFit.nNoiseDraws; end
predictionDraws=p.Results.PredictionDraws;
if isempty(predictionDraws), predictionDraws=primaryFit.nNoiseDraws; end

T=cache.trials; mask=makeMask(T,primaryFit.trialMask);
exactMask=mask & string(T.domain)==string(p.Results.Domain) & ...
  string(T.animal)==string(p.Results.Animal) & ...
  T.probeOffsetDeg==double(p.Results.ProbeOffsetDeg) & ...
  T.stepMagnitudePC==double(p.Results.StepMagnitudePC);
idx=find(exactMask); nSessions=numel(unique(string(T.sessionID(idx))));
assert(numel(idx)>=p.Results.MinimumTrials, ...
  'Only %d exact-condition trials; at least %d are required.', ...
  numel(idx),p.Results.MinimumTrials);
assert(nSessions>=p.Results.MinimumSessions, ...
  'Only %d exact-condition sessions; at least %d are required.', ...
  nSessions,p.Results.MinimumSessions);
assert(all(isfinite(primaryResult.prediction.pRDK1(idx))), ...
  'Primary prediction does not cover every selected trial.');

cfg=primaryFit.config;
cfg=maxpoolfit.configureDecisionRule(cfg,"directional-change-max");
cfg.adaptation=0; cfg.adaptationSelection="directional-change optimum at boundary";
cfg.verbose=false; cfg.nNoiseDraws=fitDraws; cfg.finalNoiseDraws=predictionDraws;
directionalFit=maxpoolfit.fitModel(cache,cfg,'TrialMask',mask,'NDraws',fitDraws);
directionalPrediction=maxpoolfit.predict(cache,directionalFit,cfg, ...
  'TrialMask',exactMask,'NDraws',predictionDraws);

currentD=maxpoolfit.noiseCoefficientsByStep(cache,primaryResult.prediction, ...
  'TrialMask',exactMask,'MinimumTrials',p.Results.MinimumTrials, ...
  'ReliableSessionCount',p.Results.MinimumSessions);
directionalD=maxpoolfit.noiseCoefficientsByStep(cache,directionalPrediction, ...
  'TrialMask',exactMask,'MinimumTrials',p.Results.MinimumTrials, ...
  'ReliableSessionCount',p.Results.MinimumSessions);
observed=oneRow(currentD.summary,"observed");
current=oneRow(currentD.summary,"primary domain noise");
directional=oneRow(directionalD.summary,"primary domain noise");
observed.source="observed";
current.source="absolute-step-max";
directional.source="directional-change-max";
gainSummary=[observed;current;directional];

model=["absolute-step-max";"directional-change-max"];
adaptation=[primaryFit.adaptation;directionalFit.adaptation];
logLikelihood=[primaryFit.logLikelihood;directionalFit.logLikelihood];
deltaLogLFromAbsolute=logLikelihood-logLikelihood(1);
nTrials=[primaryFit.nTrials;directionalFit.nTrials];
nNoiseDraws=[primaryFit.nNoiseDraws;fitDraws];
modelComparison=table(model,adaptation,nTrials,nNoiseDraws, ...
  logLikelihood,deltaLogLFromAbsolute);
currentDomain=primaryFit.domainFit; currentDomain.model= ...
  repmat("absolute-step-max",height(currentDomain),1);
directionalDomain=directionalFit.domainFit; directionalDomain.model= ...
  repmat("directional-change-max",height(directionalDomain),1);
domainComparison=[currentDomain;directionalDomain];

selection=struct('domain',string(p.Results.Domain), ...
  'animal',string(p.Results.Animal), ...
  'probeOffsetDeg',double(p.Results.ProbeOffsetDeg), ...
  'stepMagnitudePC',double(p.Results.StepMagnitudePC), ...
  'nTrials',numel(idx),'nSessions',nSessions);
diagnostic=struct('gainSummary',gainSummary, ...
  'modelComparison',modelComparison,'domainComparison',domainComparison, ...
  'selection',selection,'trialIdx',idx, ...
  'directionalFit',directionalFit,'directionalPrediction',directionalPrediction, ...
  'fitDraws',fitDraws,'predictionDraws',predictionDraws, ...
  'createdAt',datetime('now'));
end

function row=oneRow(S,source)
row=S(string(S.source)==string(source),:);
assert(height(row)==1,'Exact-condition diagnostic did not resolve to one %s row.',source);
end

function mask=makeMask(T,input)
if islogical(input)
  assert(numel(input)==height(T),'TrialMask must contain one value per cached trial.');
  mask=input(:);
else
  mask=false(height(T),1); mask(input)=true;
end
end
