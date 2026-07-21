function poolingFits = fitIDQOpponentPooling(trialTable, sessionFits, ...
  alignedWeibull, targetPerformance, preferredToNullRatio)
% fitIDQOpponentPooling
%
% Compare summation and max pooling after a fixed opponent rectification.
% For each signed stream n, the two positive-going mechanism increments are
%   preferred: max(n,0)
%   opponent:  max(-n,0)/preferredToNullRatio.
%
% The primary comparison leaves the change side empirical:
%   cEff = cStep + gD*nChangeD + gCN*(nChangePlus+nChangeMinus)
%          + gPool*poolNoChange.
% Sum and max models therefore have the same three fitted parameters and
% differ only in no-change-side pooling.
%
% A secondary 2-by-2 comparison applies sum or max pooling independently
% to both sides, with separate unconstrained change- and no-change-side
% gains. Pool predictors are centered within session. Positive no-change
% activation retains a positive predictor sign, so its fitted gain is
% expected to be negative.

if nargin < 5 || isempty(preferredToNullRatio)
  preferredToNullRatio = 2.42;
end
if ~isscalar(preferredToNullRatio) || ~isfinite(preferredToNullRatio) || ...
    preferredToNullRatio <= 0
  error('fitIDQOpponentPooling:BadRatio', ...
    'preferredToNullRatio must be a finite positive scalar.');
end

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
  "dirIndex", "changeNoisePredDir1", "changeNoisePredDir2", ...
  "changeNoisePredDir3", "noChangeNoisePredDir1", ...
  "noChangeNoisePredDir2", "noChangeNoisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
  error('fitIDQOpponentPooling:MissingVariables', ...
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
  error('fitIDQOpponentPooling:MissingSessionThreshold', ...
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
nNoChangeAll = sum(XnoChange, 2);
aNoChange = centerWithinSession(sum(abs(XnoChange), 2), sessionIndex);

% Six candidate activations per side: three preferred and three opponent.
[changeSumRaw, changeMaxRaw, changeWinner] = makePools( ...
  Xchange, preferredToNullRatio, stepCoh, dirIndex);
[noChangeSumRaw, noChangeMaxRaw, noChangeWinner] = makePools( ...
  XnoChange, preferredToNullRatio, zeros(nTrials, 1), dirIndex);

% Remove the deterministic step from change-side activation before
% centering. Centering keeps these tests focused on trialwise correspondence
% rather than a mean noise/no-noise performance shift.
changeSum = centerWithinSession(changeSumRaw - stepCoh, sessionIndex);
changeMax = centerWithinSession(changeMaxRaw - stepCoh, sessionIndex);
noChangeSum = centerWithinSession(noChangeSumRaw, sessionIndex);
noChangeMax = centerWithinSession(noChangeMaxRaw, sessionIndex);

common = struct('stepCoh', stepCoh, 'correct', correct, ...
  'sessionAlpha', sessionAlpha, 'betaWeibull', betaWeibull, 'lapse', lapse);

empiricalLinear = fitGainModel('empiricalLinear', ...
  ["g_D", "g_ChangeN", "g_NoChangeLinear"], ...
  [nChangeD, nChangeN, nNoChangeAll], [1.3; -0.1; -0.1], common);
empiricalFreeSigned = fitGainModel('empiricalFreeSigned', ...
  ["g_D", "g_ChangeN", "g_NoChange", "h_NoChange"], ...
  [nChangeD, nChangeN, nNoChangeAll, aNoChange], ...
  [empiricalLinear.gain; -0.3], common);
controlledSum = fitGainModel('controlledNoChangeSum', ...
  ["g_D", "g_ChangeN", "g_NoChangePool"], ...
  [nChangeD, nChangeN, noChangeSum], [1.3; -0.1; -0.4], common);
controlledMax = fitGainModel('controlledNoChangeMax', ...
  ["g_D", "g_ChangeN", "g_NoChangePool"], ...
  [nChangeD, nChangeN, noChangeMax], controlledSum.gain, common);

poolingNames = {'changeSum_noChangeSum', 'changeSum_noChangeMax', ...
  'changeMax_noChangeSum', 'changeMax_noChangeMax'};
changePredictors = {changeSum, changeSum, changeMax, changeMax};
noChangePredictors = {noChangeSum, noChangeMax, noChangeSum, noChangeMax};
pooling2x2 = struct();
for iModel = 1:numel(poolingNames)
  name = poolingNames{iModel};
  pooling2x2.(name) = fitGainModel(name, ...
    ["g_ChangePool", "g_NoChangePool"], ...
    [changePredictors{iModel}, noChangePredictors{iModel}], ...
    [1; -0.4], common);
end

controlledModels = {empiricalLinear, empiricalFreeSigned, ...
  controlledSum, controlledMax};
controlledModels = addInformationCriteria(controlledModels);
empiricalLinear = controlledModels{1};
empiricalFreeSigned = controlledModels{2};
controlledSum = controlledModels{3};
controlledMax = controlledModels{4};

twoByTwoModels = cell(size(poolingNames));
for iModel = 1:numel(poolingNames)
  twoByTwoModels{iModel} = pooling2x2.(poolingNames{iModel});
end
twoByTwoModels = addInformationCriteria(twoByTwoModels);
for iModel = 1:numel(poolingNames)
  pooling2x2.(poolingNames{iModel}) = twoByTwoModels{iModel};
end

% Equal-parameter, non-nested model contrasts. Positive deltaNLL favors max.
controlledContrast = makeNLLContrast(controlledSum, controlledMax, ...
  sessionIndex, correct);
changePoolingContrast = makeNLLContrast( ...
  pooling2x2.changeSum_noChangeMax, pooling2x2.changeMax_noChangeMax, ...
  sessionIndex, correct);
noChangePoolingContrast = makeNLLContrast( ...
  pooling2x2.changeMax_noChangeSum, pooling2x2.changeMax_noChangeMax, ...
  sessionIndex, correct);

winnerStats = makeWinnerStats(changeWinner, noChangeWinner, dirIndex);
poolCorrelation = corr([changeSum, changeMax, noChangeSum, noChangeMax], ...
  'Rows', 'complete');

poolingFits = struct();
poolingFits.preferredToNullRatio = preferredToNullRatio;
poolingFits.empiricalLinear = empiricalLinear;
poolingFits.empiricalFreeSigned = empiricalFreeSigned;
poolingFits.controlledSum = controlledSum;
poolingFits.controlledMax = controlledMax;
poolingFits.controlledContrast = controlledContrast;
poolingFits.pooling2x2 = pooling2x2;
poolingFits.pooling2x2Order = poolingNames;
poolingFits.changePoolingContrast = changePoolingContrast;
poolingFits.noChangePoolingContrast = noChangePoolingContrast;
poolingFits.winnerStats = winnerStats;
poolingFits.poolCorrelation = poolCorrelation;
poolingFits.nTrials = nTrials;
poolingFits.betaWeibull = betaWeibull;
poolingFits.lapse = lapse;
poolingFits.targetPerformance = targetPerformance;
poolingFits.trialData = table(sessionIndex, stepCoh, correct, dirIndex, ...
  nChangeD, nChangePlus, nChangeMinus, nChangeN, nNoChangeAll, ...
  changeSumRaw, changeMaxRaw, noChangeSumRaw, noChangeMaxRaw, ...
  changeSum, changeMax, noChangeSum, noChangeMax, ...
  changeWinner, noChangeWinner, empiricalLinear.pPred, ...
  empiricalFreeSigned.pPred, controlledSum.pPred, controlledMax.pPred, ...
  'VariableNames', {'sessionIndex', 'stepCoh', 'correct', 'dirIndex', ...
  'nChangeD', 'nChangePlus', 'nChangeMinus', 'nChangeN', 'nNoChangeAll', ...
  'changeSumRaw', 'changeMaxRaw', 'noChangeSumRaw', 'noChangeMaxRaw', ...
  'changeSum', 'changeMax', 'noChangeSum', 'noChangeMax', ...
  'changeWinner', 'noChangeWinner', 'pEmpiricalLinear', ...
  'pEmpiricalFreeSigned', 'pControlledSum', 'pControlledMax'});
poolingFits.predictorDefinition = sprintf([ ...
  'Preferred=max(n,0); opponent=max(-n,0)/%.6g; sum or max across ' ...
  'six candidates; nonlinear pool predictors centered within session.'], ...
  preferredToNullRatio);

end

%% ------------------------------------------------------------------------
function [sumPool, maxPool, winner] = makePools(X, ratio, stepCoh, dirIndex)

XwithStep = X;
row = (1:size(X, 1))';
linearIndex = sub2ind(size(X), row, dirIndex);
XwithStep(linearIndex) = XwithStep(linearIndex) + stepCoh;
preferred = max(XwithStep, 0);
opponent = max(-XwithStep, 0) / ratio;
candidates = [preferred, opponent];
sumPool = sum(candidates, 2);
[maxPool, winner] = max(candidates, [], 2);

end

%% ------------------------------------------------------------------------
function stats = makeWinnerStats(changeWinner, noChangeWinner, dirIndex)

nTrials = numel(changeWinner);
row = (1:nTrials)';
changeWinnerDirection = mod(changeWinner - 1, 3) + 1;
noChangeWinnerDirection = mod(noChangeWinner - 1, 3) + 1;
stats = struct();
stats.changePreferredFraction = mean(changeWinner <= 3);
stats.changeOpponentFraction = mean(changeWinner > 3);
stats.changeDriftMatchedFraction = mean(changeWinnerDirection == dirIndex);
stats.changeNonDriftFraction = 1 - stats.changeDriftMatchedFraction;
stats.noChangePreferredFraction = mean(noChangeWinner <= 3);
stats.noChangeOpponentFraction = mean(noChangeWinner > 3);
stats.noChangeRelativeDirectionFraction = zeros(1, 3);
relativeDirection = mod(noChangeWinnerDirection - dirIndex, 3) + 1;
for iDir = 1:3
  stats.noChangeRelativeDirectionFraction(iDir) = mean(relativeDirection == iDir);
end
stats.changeWinner = changeWinner;
stats.noChangeWinner = noChangeWinner;
stats.row = row;

end

%% ------------------------------------------------------------------------
function contrast = makeNLLContrast(sumFit, maxFit, sessionIndex, correct)

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
contrast.nSessionsFavoringSum = sum(sessionDeltaNLL < 0);

end

%% ------------------------------------------------------------------------
function p = clampProbability(p)
p = min(max(p, eps), 1 - eps);
end

%% ------------------------------------------------------------------------
function models = addInformationCriteria(models)

for iModel = 1:numel(models)
  models{iModel}.AIC = 2 * models{iModel}.nParameters + ...
    2 * models{iModel}.negLogLikelihood;
end
bestAIC = min(cellfun(@(x) x.AIC, models));
for iModel = 1:numel(models)
  models{iModel}.deltaAIC = models{iModel}.AIC - bestAIC;
end

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
function fit = fitGainModel(modelName, predictorNames, X, gain0, common)

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

hessian = finiteDifferenceHessian(objective, gainHat);
SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);
covariance = nan(nParameters);
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
