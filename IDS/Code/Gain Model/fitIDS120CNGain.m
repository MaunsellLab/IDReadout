function gainFits = fitIDS120CNGain(trialTable, sessionFits, alignedWeibull, targetPerformance, varargin)
% fitIDS120CNGain  Fit step-noise gains for both RDKs in 120-degree IDQ/IDS.
%
% Predictor order in the unconstrained fit is:
%   C drift, C +probe, C -probe, N drift, N +probe, N -probe.
%
% The supplied predictors must already be weighted by a unit-integral,
% leave-one-session-out step kernel.  Prestep noise is not fitted.

p = inputParser;
p.addParameter('NBootstrap',0,@(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('BootstrapSeed',3101,@(x) isnumeric(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

required = {'sessionIndex','stepCoh','correct','hasStepNoise', ...
  'noisePredCDrift','noisePredCPlus','noisePredCMinus', ...
  'noisePredNDrift','noisePredNPlus','noisePredNMinus'};
missing = required(~ismember(required, trialTable.Properties.VariableNames));
if ~isempty(missing)
  error('fitIDS120CNGain:MissingVariables', ...
    'Trial table lacks required variable(s): %s.', strjoin(missing, ', '));
end

idx = logical(trialTable.hasStepNoise);
sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
X6 = double([trialTable.noisePredCDrift(idx), trialTable.noisePredCPlus(idx), ...
  trialTable.noisePredCMinus(idx), trialTable.noisePredNDrift(idx), ...
  trialTable.noisePredNPlus(idx), trialTable.noisePredNMinus(idx)]);
valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & ...
  all(isfinite(X6), 2);
sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
X6 = X6(valid,:);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;
sessionThreshold = nan(size(sessionIndex));
for k = 1:height(sessionFits)
  sessionThreshold(sessionIndex == sessionFits.sessionIndex(k)) = sessionFits.threshold(k);
end
if any(~isfinite(sessionThreshold))
  error('fitIDS120CNGain:MissingSessionThreshold', ...
    'Could not assign a threshold to every noisy trial.');
end
sessionAlpha = idqWeibullAlphaForThreshold( ...
  targetPerformance, sessionThreshold, betaWeibull, lapse);

common = struct('stepCoh', stepCoh, 'correct', correct, ...
  'sessionIndex', sessionIndex, 'sessionAlpha', sessionAlpha, ...
  'betaWeibull', betaWeibull, 'lapse', lapse);

XSideRole = [X6(:,1), X6(:,2)+X6(:,3), X6(:,4), X6(:,5)+X6(:,6)];
XCNP = [X6(:,1), X6(:,2)+X6(:,3), sum(X6(:,4:6),2)];
XSide = [sum(X6(:,1:3),2), sum(X6(:,4:6),2)];
[nProbeInteraction,interactionScaling] = standardizedInteraction( ...
  X6(:,5),X6(:,6),sessionIndex);
XSideRoleInteraction = [XSideRole,nProbeInteraction];

gainFits = struct();
gainFits.full = fitGainModel('full', ...
  ["C Drift","C +Probe","C -Probe","N Drift","N +Probe","N -Probe"], ...
  X6, [1; 0; 0; 0; 0; 0], common);
gainFits.sideRole = fitGainModel('sideRole', ...
  ["C Drift","C Probe","N Drift","N Probe"], ...
  XSideRole, [1; 0; 0; 0], common);
gainFits.noChangeProbeInteraction = fitGainModel('noChangeProbeInteraction', ...
  ["C Drift","C Probe","N Drift","N Probe","N Probe Interaction"], ...
  XSideRoleInteraction,[gainFits.sideRole.gain;0],common);
gainFits.noChangeProbeInteraction.interactionScaling = interactionScaling;
gainFits.noChangeProbeInteraction.comparison = nestedComparison( ...
  gainFits.sideRole,gainFits.noChangeProbeInteraction,1);
gainFits.cnp = fitGainModel('cnp', ["C Drift","C Probe","N All"], ...
  XCNP, [1; 0; 0], common);
gainFits.side = fitGainModel('side', ["C All","N All"], ...
  XSide, [1; 0], common);
gainFits.sideRole = addSessionBootstrap(gainFits.sideRole,XSideRole,common, ...
  round(opts.NBootstrap),opts.BootstrapSeed);
gainFits.noChangeProbeInteraction = addSessionBootstrap( ...
  gainFits.noChangeProbeInteraction,XSideRoleInteraction,common, ...
  round(opts.NBootstrap),opts.BootstrapSeed+10000);
gainFits.symmetry = makeSymmetryTests(gainFits.full);
gainFits.contrasts = struct();
gainFits.contrasts.changeProbeMinusNoChangeProbe = gainContrast( ...
  gainFits.sideRole,[0 1 0 -1],0,'C Probe - N Probe');
gainFits.predictorOrder = gainFits.full.predictorNames;
end

function [interaction,info] = standardizedInteraction(xPlus,xMinus,sessionIndex)
% A local test for nonadditivity between the two noChange probes.  Center
% each predictor within session, scale each centered predictor by its
% pooled SD, multiply, and center the product within session.  The fitted
% coefficient is therefore an effective-coherence change for a joint
% one-SD modulation.  A negative value is the signature expected from
% mutual suppression around the operating point.
xPlusCentered = centerWithinSession(double(xPlus),sessionIndex);
xMinusCentered = centerWithinSession(double(xMinus),sessionIndex);
scalePlus = std(xPlusCentered,0,'omitnan');
scaleMinus = std(xMinusCentered,0,'omitnan');
if ~isfinite(scalePlus) || ~isfinite(scaleMinus) || ...
    scalePlus <= 0 || scaleMinus <= 0
  error('fitIDS120CNGain:BadInteractionScaling', ...
    'NoChange probe predictors have invalid within-session SDs.');
end
zPlus = xPlusCentered/scalePlus;
zMinus = xMinusCentered/scaleMinus;
interaction = centerWithinSession(zPlus.*zMinus,sessionIndex);
info = struct('definition', ...
  'within-session centered product of standardized N +Probe and N -Probe predictors', ...
  'plusPredictorSD',scalePlus,'minusPredictorSD',scaleMinus, ...
  'centeredPredictorCorrelation',corr(zPlus,zMinus,'Rows','complete'), ...
  'interactionMean',mean(interaction,'omitnan'), ...
  'interactionSD',std(interaction,0,'omitnan'));
end

function y = centerWithinSession(x,sessionIndex)
y = nan(size(x));
sessions = unique(sessionIndex);
for k = 1:numel(sessions)
  use = sessionIndex == sessions(k);
  y(use) = x(use)-mean(x(use),'omitnan');
end
end

function C = nestedComparison(nullFit,alternativeFit,df)
deltaNLL = nullFit.negLogLikelihood-alternativeFit.negLogLikelihood;
statistic = max(0,2*deltaNLL);
C = struct('nullModel',nullFit.model,'alternativeModel',alternativeFit.model, ...
  'deltaNegLogLikelihood',deltaNLL,'likelihoodRatio',statistic,'df',df, ...
  'pValue',1-gammainc(statistic/2,df/2));
end

function fit = fitGainModel(modelName, predictorNames, X, gain0, common)
nParameters = size(X,2);
lb = -5 * ones(nParameters,1);
ub = 5 * ones(nParameters,1);
objective = @(gain) gainNLL(gain, X, common);
opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
  'OptimalityTolerance', 1e-8, 'StepTolerance', 1e-10, ...
  'MaxFunctionEvaluations', 10000, 'MaxIterations', 3000);
[gainHat,nll,exitflag] = fmincon(objective,gain0(:),[],[],[],[],lb,ub,[],opts);
gainHat = gainHat(:);

hessian = finiteDifferenceHessian(objective,gainHat);
[hessianCovariance, hessianOK] = covarianceFromHessian(hessian);
clusterCovariance = clusterRobustCovariance(gainHat,X,common,hessianCovariance,hessianOK);
if all(isfinite(clusterCovariance),'all')
  covariance = clusterCovariance;
  covarianceMethod = 'Session-cluster robust sandwich covariance';
else
  covariance = hessianCovariance;
  covarianceMethod = 'Observed-Hessian covariance';
end
[SE,CI95] = covarianceIntervals(gainHat,covariance);

fit = struct('model',modelName,'predictorNames',predictorNames(:)', ...
  'gain',gainHat,'SE',SE,'CI95',CI95,'covariance',covariance, ...
  'covarianceMethod',covarianceMethod,'hessianCovariance',hessianCovariance, ...
  'clusterCovariance',clusterCovariance,'hessian',hessian, ...
  'negLogLikelihood',nll,'exitflag',exitflag,'nParameters',nParameters, ...
  'nTrials',numel(common.correct),'nSessions',numel(unique(common.sessionIndex)), ...
  'nEffectiveCohClipped',sum(common.stepCoh + X*gainHat < 0));
end

function nll = gainNLL(gain, X, common)
nll = weightedGainNLL(gain,X,common,ones(size(common.correct)));
end

function nll = weightedGainNLL(gain,X,common,trialWeight)
effectiveCoh = max(common.stepCoh + X*gain(:),0);
p = idqWeibullP(effectiveCoh,common.sessionAlpha,common.betaWeibull,common.lapse);
p = min(max(p,eps),1-eps);
nll = -sum(trialWeight.*(common.correct.*log(p) + ...
  (1-common.correct).*log(1-p)));
end

function fit = addSessionBootstrap(fit,X,common,nBootstrap,seed)
fit.bootstrapNRequested = nBootstrap;
fit.bootstrapNSuccess = 0;
fit.bootstrapGain = nan(nBootstrap,numel(fit.gain));
fit.bootstrapCI95 = nan(numel(fit.gain),2);
if nBootstrap <= 0
  return
end

sessions = unique(common.sessionIndex);
nSessions = numel(sessions);
if nSessions < 2
  return
end
[~,sessionLocation] = ismember(common.sessionIndex,sessions);
oldState = rng;
cleanup = onCleanup(@() rng(oldState)); %#ok<NASGU>
rng(seed,'twister');
options = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-7,'StepTolerance',1e-9, ...
  'MaxFunctionEvaluations',4000,'MaxIterations',1000);
lb = -5*ones(numel(fit.gain),1);
ub = 5*ones(numel(fit.gain),1);
for b = 1:nBootstrap
  sampledLocations = randi(nSessions,nSessions,1);
  sessionCounts = accumarray(sampledLocations,1,[nSessions 1]);
  trialWeight = sessionCounts(sessionLocation);
  objective = @(gain) weightedGainNLL(gain,X,common,trialWeight);
  [gainHat,~,exitflag] = fmincon(objective,fit.gain,[],[],[],[], ...
    lb,ub,[],options);
  if exitflag > 0 && all(isfinite(gainHat))
    fit.bootstrapGain(b,:) = gainHat(:)';
  end
end
valid = all(isfinite(fit.bootstrapGain),2);
fit.bootstrapNSuccess = sum(valid);
for k = 1:numel(fit.gain)
  fit.bootstrapCI95(k,:) = empiricalCI(fit.bootstrapGain(valid,k));
end
end

function C = gainContrast(fit,R,r,label)
R = R(:)';
estimate = R*fit.gain-r;
variance = R*fit.covariance*R';
SE = nan;
CI95 = [nan nan];
pValue = nan;
if isfinite(variance) && variance > 0
  SE = sqrt(variance);
  CI95 = estimate + [-1 1].*(1.96*SE);
  pValue = erfc(abs(estimate/SE)/sqrt(2));
end
bootstrapValues = [];
bootstrapCI95 = [nan nan];
bootstrapPValue = nan;
if isfield(fit,'bootstrapGain') && ~isempty(fit.bootstrapGain)
  valid = all(isfinite(fit.bootstrapGain),2);
  bootstrapValues = fit.bootstrapGain(valid,:)*R'-r;
  bootstrapCI95 = empiricalCI(bootstrapValues);
  if ~isempty(bootstrapValues)
    n = numel(bootstrapValues);
    centeredDeviation = bootstrapValues-estimate;
    bootstrapPValue = (sum(abs(centeredDeviation) >= abs(estimate))+1)/(n+1);
  end
end
C = struct('label',label,'contrast',R,'nullValue',r, ...
  'estimate',estimate,'SE',SE,'CI95',CI95,'pValue',pValue, ...
  'bootstrapValues',bootstrapValues,'bootstrapCI95',bootstrapCI95, ...
  'bootstrapPValue',bootstrapPValue);
end

function ci = empiricalCI(values)
values = sort(values(isfinite(values)));
if isempty(values)
  ci = [nan nan];
  return
end
n = numel(values);
indices = max(1,min(n,round([0.025 0.975]*(n-1)+1)));
ci = values(indices)';
end

function C = clusterRobustCovariance(gain,X,common,hessianCovariance,hessianOK)
nParameters = numel(gain);
C = nan(nParameters);
sessions = unique(common.sessionIndex);
nClusters = numel(sessions);
nTrials = numel(common.correct);
if ~hessianOK || nClusters < 2 || nTrials <= nParameters
  return
end
meat = zeros(nParameters);
for k = 1:nClusters
  use = common.sessionIndex == sessions(k);
  local = subsetCommon(common,use);
  objective = @(g) gainNLL(g,X(use,:),local);
  score = finiteDifferenceGradient(objective,gain);
  meat = meat + score*score';
end
smallSample = (nClusters/(nClusters-1))*((nTrials-1)/(nTrials-nParameters));
C = smallSample * hessianCovariance * meat * hessianCovariance;
C = (C+C')/2;
if any(~isfinite(C),'all') || any(diag(C) <= 0)
  C(:) = nan;
end
end

function out = subsetCommon(in,use)
out = in;
fields = {'stepCoh','correct','sessionIndex','sessionAlpha'};
for k = 1:numel(fields)
  out.(fields{k}) = in.(fields{k})(use);
end
end

function tests = makeSymmetryTests(fit)
tests = struct();
tests.changeProbe = waldTest(fit,[0 1 -1 0 0 0],0,'C +Probe = C -Probe');
tests.noChangeProbe = waldTest(fit,[0 0 0 0 1 -1],0,'N +Probe = N -Probe');
tests.noChangeAll = waldTest(fit,[0 0 0 1 -1 0; 0 0 0 1 0 -1], ...
  [0;0],'N Drift = N +Probe = N -Probe');
end

function T = waldTest(fit,R,r,label)
delta = R*fit.gain-r;
V = R*fit.covariance*R';
df = size(R,1);
statistic = nan;
pValue = nan;
if all(isfinite(V),'all') && rcond(V) > 1e-12
  statistic = delta'*(V\delta);
  pValue = 1-gammainc(statistic/2,df/2);
end
T = struct('label',label,'contrast',R,'estimate',delta, ...
  'statistic',statistic,'df',df,'pValue',pValue);
end

function [C,ok] = covarianceFromHessian(H)
C = nan(size(H));
ok = false;
if all(isfinite(H),'all')
  H = (H+H')/2;
  if all(eig(H) > 0)
    C = inv(H);
    ok = all(isfinite(C),'all') && all(diag(C) > 0);
  end
end
end

function [SE,CI95] = covarianceIntervals(estimate,C)
SE = nan(size(estimate));
CI95 = nan(numel(estimate),2);
if any(~isfinite(C),'all')
  return
end
variance = diag(C);
ok = isfinite(variance) & variance > 0;
SE(ok) = sqrt(variance(ok));
CI95(ok,:) = estimate(ok) + [-1 1].*(1.96*SE(ok));
end

function gradient = finiteDifferenceGradient(fun,x)
x = x(:);
h = 3e-5.*max(1,abs(x));
gradient = nan(size(x));
for k = 1:numel(x)
  e = zeros(size(x));
  e(k) = h(k);
  gradient(k) = (fun(x+e)-fun(x-e))/(2*h(k));
end
end

function H = finiteDifferenceHessian(fun,x)
x = x(:);
n = numel(x);
H = nan(n);
h = 3e-4.*max(1,abs(x));
f0 = fun(x);
for i = 1:n
  ei = zeros(n,1); ei(i) = h(i);
  H(i,i) = (fun(x+ei)-2*f0+fun(x-ei))/h(i)^2;
  for j = i+1:n
    ej = zeros(n,1); ej(j) = h(j);
    H(i,j) = (fun(x+ei+ej)-fun(x+ei-ej)-fun(x-ei+ej)+fun(x-ei-ej)) / ...
      (4*h(i)*h(j));
    H(j,i) = H(i,j);
  end
end
end
