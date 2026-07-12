function gainFits = fitIDQNoiseGain(trialTable, sessionFits, alignedWeibull, targetPerformance)
% fitIDQNoiseGain
%
% Fit shared gains for four IDQ directional-noise models:
%   combined         one gain for the sum of all three streams
%   driftRelative    separate gains for drift, +120, and -120 streams
%   absolute         separate gains for absolute direction streams 1:3
%   driftNonDrift    separate gains for drift and pooled non-drift streams
%
% Required trialTable variables:
%   sessionIndex, stepCoh, correct, hasStepNoise, dirIndex,
%   noisePredDir1, noisePredDir2, noisePredDir3

requiredVars = ["sessionIndex", "stepCoh", "correct", "hasStepNoise", ...
    "dirIndex", "noisePredDir1", "noisePredDir2", "noisePredDir3"];
missingVars = setdiff(requiredVars, string(trialTable.Properties.VariableNames));
if ~isempty(missingVars)
    error('fitIDQNoiseGain:MissingVariables', ...
        'trialTable is missing required variables: %s', strjoin(missingVars, ', '));
end

idx = trialTable.hasStepNoise;

sessionIndex = double(trialTable.sessionIndex(idx));
stepCoh = double(trialTable.stepCoh(idx));
correct = double(trialTable.correct(idx));
dirIndex = double(trialTable.dirIndex(idx));
Xabs = double([trialTable.noisePredDir1(idx), ...
               trialTable.noisePredDir2(idx), ...
               trialTable.noisePredDir3(idx)]);

valid = isfinite(sessionIndex) & isfinite(stepCoh) & isfinite(correct) & ...
    isfinite(dirIndex) & all(isfinite(Xabs), 2) & ...
    dirIndex >= 1 & dirIndex <= 3 & dirIndex == round(dirIndex);

sessionIndex = sessionIndex(valid);
stepCoh = stepCoh(valid);
correct = correct(valid);
dirIndex = dirIndex(valid);
Xabs = Xabs(valid, :);

betaWeibull = alignedWeibull.betaWeibull;
lapse = alignedWeibull.lapse;

sessionThreshold = nan(size(sessionIndex));
for i = 1:height(sessionFits)
    idxSession = sessionIndex == sessionFits.sessionIndex(i);
    sessionThreshold(idxSession) = sessionFits.threshold(i);
end

if any(~isfinite(sessionThreshold))
    error('fitIDQNoiseGain:MissingSessionThreshold', ...
        'Could not assign threshold to all noisy trials.');
end

sessionAlpha = idqWeibullAlphaForThreshold( ...
    targetPerformance, sessionThreshold, betaWeibull, lapse);

nTrials = numel(correct);
row = (1:nTrials)';
plusDir = mod(dirIndex, 3) + 1;
minusDir = mod(dirIndex - 2, 3) + 1;

xDrift = Xabs(sub2ind(size(Xabs), row, dirIndex));
xPlus120 = Xabs(sub2ind(size(Xabs), row, plusDir));
xMinus120 = Xabs(sub2ind(size(Xabs), row, minusDir));

Xcombined = sum(Xabs, 2);
Xrelative = [xDrift, xPlus120, xMinus120];
XdriftNonDrift = [xDrift, xPlus120 + xMinus120];

common = struct();
common.stepCoh = stepCoh;
common.correct = correct;
common.sessionAlpha = sessionAlpha;
common.betaWeibull = betaWeibull;
common.lapse = lapse;

gainFits = struct();
gainFits.combined = fitGainModel( ...
    'combined', ["allDirections"], Xcombined, 1/3, common);
gainFits.driftRelative = fitGainModel( ...
    'driftRelative', ["drift", "plus120", "minus120"], ...
    Xrelative, [1; 0; 0], common);
gainFits.absolute = fitGainModel( ...
    'absolute', ["dir1", "dir2", "dir3"], ...
    Xabs, repmat(1/3, 3, 1), common);
gainFits.driftNonDrift = fitGainModel( ...
    'driftNonDrift', ["drift", "nonDrift"], ...
    XdriftNonDrift, [1; 0], common);

% 
% 
% fprintf('SD drift:      %.4f\n', std(xDrift));
% fprintf('SD +120:       %.4f\n', std(xPlus120));
% fprintf('SD -120:       %.4f\n', std(xMinus120));
% fprintf('SD pooled non: %.4f\n', std(xPlus120 + xMinus120));
% fprintf('Corr non-drift: %.4f\n', corr(xPlus120, xMinus120));

end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName, predictorNames, X, gain0, common)

nParameters = size(X, 2);
lb = -5 * ones(nParameters, 1);
ub =  5 * ones(nParameters, 1);
gain0 = gain0(:);

objective = @(gain) negLogLikelihoodGain(gain, X, common.stepCoh, ...
    common.correct, common.sessionAlpha, common.betaWeibull, common.lapse);

opts = optimoptions('fmincon', ...
    'Display', 'off', ...
    'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 5000, ...
    'MaxIterations', 2000);

% H = hessian;       % or temporarily capture it
% disp(modelName);
% disp(predictorNames);
% C = inv(H);
% H
% C
% eig(H)
% cond(H)
% sqrt(diag(C))

% [gainHat, nll, exitflag, ~, ~, ~, hessian] = ...
%     fmincon(objective, gain0, [], [], [], [], lb, ub, [], opts);
% gainHat = gainHat(:);
% SE = nan(nParameters, 1);
% CI95 = nan(nParameters, 2);
% 
% if all(isfinite(hessian), 'all')
%     covariance = pinv(hessian);
%     variance = diag(covariance);
%     validVariance = isfinite(variance) & variance > 0;
%     SE(validVariance) = sqrt(variance(validVariance));
%     CI95(validVariance, :) = gainHat(validVariance) + ...
%         [-1 1] .* (1.96 * SE(validVariance));
% end

[gainHat, nll, exitflag] = ...
  fmincon(objective, gain0, [], [], [], [], lb, ub, [], opts);

gainHat = gainHat(:);

% Compute the observed Hessian directly from the objective at the optimum.
hessian = finiteDifferenceHessian(objective, gainHat);

SE = nan(nParameters, 1);
CI95 = nan(nParameters, 2);

if all(isfinite(hessian), 'all')
  hessian = (hessian + hessian') / 2;

  eigenvalues = eig(hessian);

  if all(eigenvalues > 0)
    covariance = inv(hessian);
    variance = diag(covariance);

    validVariance = isfinite(variance) & variance > 0;
    SE(validVariance) = sqrt(variance(validVariance));

    CI95(validVariance, :) = gainHat(validVariance) + ...
      [-1 1] .* (1.96 * SE(validVariance));
  end
end

fit = struct();
fit.model = modelName;
fit.predictorNames = predictorNames(:)';
fit.gain = gainHat;
fit.SE = SE;
fit.CI95 = CI95;
fit.negLogLikelihood = nll;
fit.exitflag = exitflag;
fit.nParameters = nParameters;
fit.nTrials = numel(common.correct);
fit.nEffectiveCohClipped = sum(common.stepCoh + X * gainHat < 0);

% the following can be removed after diagnostics are settled.

% fit.hessian = hessian;
% fit.hessianEigenvalues = eig(hessian);
% fit.hessianCondition = cond(hessian);
% 
% disp(fit.model);
% disp(fit.gain);
% disp(fit.SE);
% disp(fit.hessian);
% disp(eig(fit.hessian));
% disp(cond(fit.hessian));

end

%% ------------------------------------------------------------------------
function nll = negLogLikelihoodGain(gain, X, stepCoh, correct, ...
    sessionAlpha, betaWeibull, lapse)

effectiveCoh = stepCoh + X * gain(:);

% Guard against numerical excursions during the bounded gain search.
effectiveCoh = max(effectiveCoh, 0);

p = idqWeibullP(effectiveCoh, sessionAlpha, betaWeibull, lapse);
p = min(max(p, eps), 1 - eps);

nll = -sum(correct .* log(p) + (1 - correct) .* log(1 - p));

end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(objective, x)

x = x(:);
nParameters = numel(x);

H = nan(nParameters);

% Scale-aware finite-difference step.
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

    H(i, j) = ...
      (fpp - fpm - fmp + fmm) / ...
      (4 * h(i) * h(j));

    H(j, i) = H(i, j);
  end
end

end