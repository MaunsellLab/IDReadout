function pNormFits = fitIDQOpponentPNorm(trialTable, sessionFits, ...
  alignedWeibull, targetPerformance, preferredToNullRatio)
% fitIDQOpponentPNorm
%
% Profile opponent-rectified p-norm pooling on the no-change side:
%   P_i = max(n_i,0)
%   O_i = max(-n_i,0)/preferredToNullRatio
%   M_p = (sum_i P_i^p + sum_i O_i^p)^(1/p)
% with p=1 giving summation and p=Inf giving a hard maximum.
%
% Two model contexts are fit:
%   mechanistic: hard max on change side, p-norm on no-change side
%   controlled: empirical change drift/non-drift terms, p-norm on
%               no-change side
%
% Pool predictors are centered within session. All gains are unconstrained.

if nargin < 5 || isempty(preferredToNullRatio)
  preferredToNullRatio = 2.42;
end
if ~isscalar(preferredToNullRatio) || ~isfinite(preferredToNullRatio) || ...
    preferredToNullRatio <= 0
  error('fitIDQOpponentPNorm:BadRatio', ...
    'preferredToNullRatio must be a finite positive scalar.');
end

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "changeNoisePredDir1", "changeNoisePredDir2", ...
  "changeNoisePredDir3", "noChangeNoisePredDir1", ...
  "noChangeNoisePredDir2", "noChangeNoisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQOpponentPNorm:MissingVariables', ...
    'trialTable is missing required variables: %s', strjoin(missingVars, ', '));
end

idx = trialTable.hasStepNoise;
sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
dirIndex = double(trialTable.dirIndex(idx));
Xchange = double([trialTable.changeNoisePredDir1(idx), ...
  trialTable.changeNoisePredDir2(idx), trialTable.changeNoisePredDir3(idx)]);
XnoChange = double([trialTable.noChangeNoisePredDir1(idx), ...
  trialTable.noChangeNoisePredDir2(idx), trialTable.noChangeNoisePredDir3(idx)]);

valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & ...
  isfinite(dirIndex) & all(isfinite(Xchange), 2) & ...
  all(isfinite(XnoChange), 2) & dirIndex >= 1 & dirIndex <= 3 & ...
  dirIndex == round(dirIndex);
sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
dirIndex = dirIndex(valid);
Xchange = Xchange(valid, :);
XnoChange = XnoChange(valid, :);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;
sessionThreshold = nan(size(sessionIndex));
for iSession = 1:height(sessionFits)
  idxSession = sessionIndex == sessionFits.sessionIndex(iSession);
  sessionThreshold(idxSession) = sessionFits.threshold(iSession);
end
if any(~isfinite(sessionThreshold))
  error('fitIDQOpponentPNorm:MissingSessionThreshold', ...
    'Could not assign threshold to all noisy trials.');
end
sessionAlpha = idqWeibullAlphaForThreshold(targetPerformance, ...
  sessionThreshold, betaWeibull, lapse);

nTrials = numel(correct);
row = (1:nTrials)';
plusDir = mod(dirIndex, 3) + 1;
minusDir = mod(dirIndex - 2, 3) + 1;
nChangeD = Xchange(sub2ind(size(Xchange), row, dirIndex));
nChangePlus = Xchange(sub2ind(size(Xchange), row, plusDir));
nChangeMinus = Xchange(sub2ind(size(Xchange), row, minusDir));
nChangeN = nChangePlus + nChangeMinus;

changeCandidates = makeCandidates(Xchange, preferredToNullRatio, ...
  stepCoh, dirIndex);
noChangeCandidates = makeCandidates(XnoChange, preferredToNullRatio, ...
  zeros(nTrials, 1), dirIndex);
changeMaxRaw = max(changeCandidates, [], 2);
changeMax = centerWithinSession(changeMaxRaw - stepCoh, sessionIndex);

common = struct('stepCoh', stepCoh, 'correct', correct, ...
  'sessionAlpha', sessionAlpha, 'betaWeibull', betaWeibull, 'lapse', lapse);
data = struct('sessionIndex', sessionIndex, 'nChangeD', nChangeD, ...
  'nChangeN', nChangeN, 'changeMax', changeMax, ...
  'noChangeCandidates', noChangeCandidates);

pGridFinite = [1 1.25 1.5 2 3 4 6 8 12 16 24 32 48 64 96 128];
pUpper = pGridFinite(end);
mechanisticProfile = fitProfile('mechanistic', pGridFinite, data, common);
controlledProfile = fitProfile('controlled', pGridFinite, data, common);

mechanisticContinuous = fitContinuousP('mechanistic', pUpper, data, common);
controlledContinuous = fitContinuousP('controlled', pUpper, data, common);

mechanisticHardMax = fitAtP('mechanistic', Inf, data, common, [], true);
controlledHardMax = fitAtP('controlled', Inf, data, common, [], true);
mechanisticSum = mechanisticProfile.fits{1};
controlledSum = controlledProfile.fits{1};

mechanisticContinuous = addContinuousInformationCriteria( ...
  mechanisticContinuous, mechanisticSum, mechanisticHardMax);
controlledContinuous = addContinuousInformationCriteria( ...
  controlledContinuous, controlledSum, controlledHardMax);
mechanisticProfile = addProfileStatistics(mechanisticProfile, ...
  mechanisticContinuous, mechanisticHardMax);
controlledProfile = addProfileStatistics(controlledProfile, ...
  controlledContinuous, controlledHardMax);

mechanisticContrast = makeEndpointContrast(mechanisticSum, mechanisticHardMax, ...
  sessionIndex, correct);
controlledContrast = makeEndpointContrast(controlledSum, controlledHardMax, ...
  sessionIndex, correct);

pNormFits = struct();
pNormFits.preferredToNullRatio = preferredToNullRatio;
pNormFits.pUpper = pUpper;
pNormFits.mechanistic = struct('sum', mechanisticSum, ...
  'continuous', mechanisticContinuous, 'hardMax', mechanisticHardMax, ...
  'profile', mechanisticProfile, 'endpointContrast', mechanisticContrast);
pNormFits.controlled = struct('sum', controlledSum, ...
  'continuous', controlledContinuous, 'hardMax', controlledHardMax, ...
  'profile', controlledProfile, 'endpointContrast', controlledContrast);
pNormFits.nTrials = nTrials;
pNormFits.betaWeibull = betaWeibull;
pNormFits.lapse = lapse;
pNormFits.targetPerformance = targetPerformance;
pNormFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nChangeD, nChangePlus, nChangeMinus, nChangeN, changeMax, ...
  mechanisticContinuous.noChangePool, controlledContinuous.noChangePool, ...
  mechanisticContinuous.pPred, controlledContinuous.pPred, ...
  mechanisticHardMax.pPred, controlledHardMax.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nChangeD', 'nChangePlus', 'nChangeMinus', 'nChangeN', 'changeMax', ...
  'mechanisticPNorm', 'controlledPNorm', 'pMechanisticPNorm', ...
  'pControlledPNorm', 'pMechanisticMax', 'pControlledMax'});
pNormFits.predictorDefinition = sprintf([ ...
  'M_p=(sum(candidate^p))^(1/p), candidates=[max(n,0), ' ...
  'max(-n,0)/%.6g], centered within session; p=Inf is exact max.'], ...
  preferredToNullRatio);

end

%% ------------------------------------------------------------------------
function profile = fitProfile(context, pGrid, data, common)

fits = cell(size(pGrid));
gain0 = [];
for iP = 1:numel(pGrid)
  fits{iP} = fitAtP(context, pGrid(iP), data, common, gain0, false);
  gain0 = fits{iP}.gain;
end
profile = struct();
profile.p = pGrid(:);
profile.fits = fits;
profile.negLogLikelihood = cellfun(@(x) x.negLogLikelihood, fits(:));
profile.gain = cell2mat(cellfun(@(x) x.gain(:)', fits(:), ...
  'UniformOutput', false));

end

%% ------------------------------------------------------------------------
function result = fitContinuousP(context, pUpper, data, common)

qBounds = [0 log(pUpper)];
objective = @(q) profileObjective(context, exp(q), data, common);
options = optimset('Display', 'off', 'TolX', 2e-4, ...
  'MaxFunEvals', 80, 'MaxIter', 80);
[qHat, ~, exitflagP] = fminbnd(objective, ...
  qBounds(1), qBounds(2), options);

% Explicitly retain a boundary solution if it is better than fminbnd's
% interior result.
candidateP = [1, exp(qHat), pUpper];
candidateFits = cell(size(candidateP));
for iCandidate = 1:numel(candidateP)
  candidateFits{iCandidate} = fitAtP(context, candidateP(iCandidate), ...
    data, common, [], false);
end
[~, iBest] = min(cellfun(@(x) x.negLogLikelihood, candidateFits));
pHat = candidateP(iBest);
result = fitAtP(context, pHat, data, common, ...
  candidateFits{iBest}.gain, true);
result.p = pHat;
result.q = log(pHat);
result.pAtLowerBound = pHat <= 1 + 1e-6;
result.pAtUpperBound = pHat >= pUpper * (1 - 1e-6);
result.exitflagP = exitflagP;
result.nParametersIncludingP = result.nParameters + 1;

end

%% ------------------------------------------------------------------------
function nll = profileObjective(context, p, data, common)

fit = fitAtP(context, p, data, common, [], false);
nll = fit.negLogLikelihood;

end

%% ------------------------------------------------------------------------
function fit = fitAtP(context, p, data, common, gain0, computeHessian)

noChangePoolRaw = pNormRows(data.noChangeCandidates, p);
noChangePool = centerWithinSession(noChangePoolRaw, data.sessionIndex);
switch context
  case 'mechanistic'
    X = [data.changeMax, noChangePool];
    names = ["g_ChangeMax", "g_NoChangePool"];
    if isempty(gain0)
      gain0 = [1.35; -0.6];
    end
  case 'controlled'
    X = [data.nChangeD, data.nChangeN, noChangePool];
    names = ["g_D", "g_ChangeN", "g_NoChangePool"];
    if isempty(gain0)
      gain0 = [1.35; -0.12; -0.6];
    end
  otherwise
    error('fitIDQOpponentPNorm:BadContext', 'Unknown context %s.', context);
end
fit = fitGainModel(sprintf('%s_pNorm', context), names, X, gain0, ...
  common, computeHessian);
fit.p = p;
fit.noChangePool = noChangePool;

end

%% ------------------------------------------------------------------------
function candidates = makeCandidates(X, ratio, stepCoh, dirIndex)

XwithStep = X;
row = (1:size(X, 1))';
linearIndex = sub2ind(size(X), row, dirIndex);
XwithStep(linearIndex) = XwithStep(linearIndex) + stepCoh;
candidates = [max(XwithStep, 0), max(-XwithStep, 0) / ratio];

end

%% ------------------------------------------------------------------------
function value = pNormRows(candidates, p)

rowMax = max(candidates, [], 2);
if isinf(p)
  value = rowMax;
  return;
end
value = zeros(size(rowMax));
idx = rowMax > 0;
scaled = candidates(idx, :) ./ rowMax(idx);
value(idx) = rowMax(idx) .* sum(scaled .^ p, 2) .^ (1 / p);

end

%% ------------------------------------------------------------------------
function continuous = addContinuousInformationCriteria(continuous, sumFit, maxFit)

continuous.AIC = 2 * continuous.nParametersIncludingP + ...
  2 * continuous.negLogLikelihood;
sumFitAIC = 2 * sumFit.nParameters + 2 * sumFit.negLogLikelihood;
maxFitAIC = 2 * maxFit.nParameters + 2 * maxFit.negLogLikelihood;
bestAIC = min([sumFitAIC, continuous.AIC, maxFitAIC]);
continuous.sumAIC = sumFitAIC;
continuous.hardMaxAIC = maxFitAIC;
continuous.deltaAIC = continuous.AIC - bestAIC;
continuous.sumDeltaAIC = sumFitAIC - bestAIC;
continuous.hardMaxDeltaAIC = maxFitAIC - bestAIC;

end

%% ------------------------------------------------------------------------
function profile = addProfileStatistics(profile, continuous, hardMax)

referenceNLL = min(continuous.negLogLikelihood, hardMax.negLogLikelihood);
profile.deltaNLL = profile.negLogLikelihood - referenceNLL;
cutoff = 0.5 * 3.84145882069412;
inside = profile.deltaNLL <= cutoff;
if any(inside)
  profile.profileCI95Low = min(profile.p(inside));
else
  profile.profileCI95Low = NaN;
end
if hardMax.negLogLikelihood - referenceNLL <= cutoff || ...
    (any(inside) && inside(end))
  profile.profileCI95High = Inf;
elseif any(inside)
  profile.profileCI95High = max(profile.p(inside));
else
  profile.profileCI95High = NaN;
end
profile.hardMaxDeltaNLL = hardMax.negLogLikelihood - referenceNLL;

end

%% ------------------------------------------------------------------------
function contrast = makeEndpointContrast(sumFit, maxFit, sessionIndex, correct)

logLikeSum = correct .* log(clampProbability(sumFit.pPred)) + ...
  (1 - correct) .* log(1 - clampProbability(sumFit.pPred));
logLikeMax = correct .* log(clampProbability(maxFit.pPred)) + ...
  (1 - correct) .* log(1 - clampProbability(maxFit.pPred));
sessions = unique(sessionIndex);
sessionDeltaNLL = nan(numel(sessions), 1);
for iSession = 1:numel(sessions)
  idx = sessionIndex == sessions(iSession);
  sessionDeltaNLL(iSession) = sum(logLikeMax(idx) - logLikeSum(idx));
end
contrast = struct();
contrast.deltaNLLMaxOverSum = sumFit.negLogLikelihood - maxFit.negLogLikelihood;
contrast.sessionIndex = sessions;
contrast.sessionDeltaNLLMaxOverSum = sessionDeltaNLL;
contrast.nSessionsFavoringMax = sum(sessionDeltaNLL > 0);

end

%% ------------------------------------------------------------------------
function xCentered = centerWithinSession(x, sessionIndex)

xCentered = x;
sessions = unique(sessionIndex);
for iSession = 1:numel(sessions)
  idx = sessionIndex == sessions(iSession);
  xCentered(idx) = x(idx) - mean(x(idx), 'omitnan');
end

end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName, predictorNames, X, gain0, common, computeHessian)

nParameters = size(X, 2);
lb = -5 * ones(nParameters, 1);
ub =  5 * ones(nParameters, 1);
objective = @(gain) negLogLikelihoodGain(gain, X, common.stepCoh, ...
  common.correct, common.sessionAlpha, common.betaWeibull, common.lapse);
opts = optimoptions('fmincon', 'Display', 'off', ...
  'Algorithm', 'interior-point', 'OptimalityTolerance', 1e-8, ...
  'StepTolerance', 1e-10, 'MaxFunctionEvaluations', 5000, ...
  'MaxIterations', 2000);
[gainHat, nll, exitflag] = fmincon(objective, gain0(:), ...
  [], [], [], [], lb, ub, [], opts);
gainHat = gainHat(:);

SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);
covariance = nan(nParameters);
hessian = nan(nParameters);
if computeHessian
  hessian = finiteDifferenceHessian(objective, gainHat);
  if all(isfinite(hessian), 'all')
    hessian = (hessian + hessian') / 2;
    if all(eig(hessian) > 0)
      covariance = inv(hessian);
      variance = diag(covariance);
      validVariance = isfinite(variance) & variance > 0;
      SE(validVariance) = sqrt(variance(validVariance));
      CI95(validVariance, :) = gainHat(validVariance) + ...
        [-1 1] .* (1.96 * SE(validVariance));
    end
  end
end

effectiveCohUnclipped = common.stepCoh + X * gainHat;
effectiveCoh = max(effectiveCohUnclipped, 0);
pPred = idqWeibullP(effectiveCoh, common.sessionAlpha, ...
  common.betaWeibull, common.lapse);

fit = struct();
fit.model = modelName;
fit.predictorNames = cellstr(string(predictorNames(:)'));
fit.gain = gainHat;
fit.SE = SE;
fit.CI95 = CI95;
fit.covariance = covariance;
fit.hessian = hessian;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.nParameters = nParameters;
fit.nTrials = numel(common.correct);
fit.nEffectiveCohClipped = sum(effectiveCohUnclipped < 0);
fit.effectiveCoh = effectiveCoh;
fit.pPred = pPred;

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihoodGain(gain, X, stepCoh, correct, ...
  sessionAlpha, betaWeibull, lapse)

effectiveCoh = max(stepCoh + X * gain(:), 0);
p = idqWeibullP(effectiveCoh, sessionAlpha, betaWeibull, lapse);
p = clampProbability(p);
nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end

%% ------------------------------------------------------------------------
function p = clampProbability(p)
p = min(max(p, eps), 1 - eps);
end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(objective, x)

x = x(:);
nParameters = numel(x);
H = nan(nParameters);
h = 3e-4 .* max(1, abs(x));
f0 = objective(x);
for i = 1:nParameters
  ei = zeros(nParameters, 1);
  ei(i) = h(i);
  fp = objective(x + ei);
  fm = objective(x - ei);
  H(i, i) = (fp - 2*f0 + fm) / h(i)^2;
  for j = i+1:nParameters
    ej = zeros(nParameters, 1);
    ej(j) = h(j);
    fpp = objective(x + ei + ej);
    fpm = objective(x + ei - ej);
    fmp = objective(x - ei + ej);
    fmm = objective(x - ei - ej);
    H(i, j) = (fpp - fpm - fmp + fmm) / (4 * h(i) * h(j));
    H(j, i) = H(i, j);
  end
end

end
