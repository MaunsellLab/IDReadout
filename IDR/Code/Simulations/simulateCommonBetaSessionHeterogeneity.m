function commonBetaSim = simulateCommonBetaSessionHeterogeneity(nSim, seed)
% simulateCommonBetaSessionHeterogeneity
% Test whether observed session-beta heterogeneity exceeds finite-sample noise.
%
% Null model:
%   Every session shares the pooled noise beta.
%   Each session retains:
%     - its actual effectiveNoisePC values
%     - its fitted pooled-model session intercept
%     - its original number of trials
%
% Simulated outcomes:
%   P(correct) = 0.5 + 0.5*logistic(alpha_session ...
%                                  + betaCommon*effectiveNoisePC)
%
% Each simulated session is then refit separately using the same model as
% fitBetaNoiseBySession. The resulting session-beta SD, random-effects tau,
% and I^2 are compared with the observed values.
%
% Usage:
%   commonBetaSim = simulateCommonBetaSessionHeterogeneity();
%   commonBetaSim = simulateCommonBetaSessionHeterogeneity(1000, 1);
%
% Reads:
%   BetaNoiseRegressionFit.mat
%   BetaNoiseSessionFits.mat
%   BetaSessionFitAnalysis.mat
%   SessionData/*.mat
%
% Saves:
%   BetaCommonSlopeSimulation.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

if nargin < 1 || isempty(nSim)
  nSim = 1000;
end
if nargin < 2 || isempty(seed)
  seed = 1;
end

validateattributes(nSim, {'numeric'}, ...
  {'scalar','integer','positive'});
validateattributes(seed, {'numeric'}, ...
  {'scalar','integer','nonnegative'});

rng(seed);

baseFolder = domainFolder(mfilename('fullpath'));
sessionFolder = fullfile(baseFolder, 'Data', 'FullSessions', 'BetaAnalysis');
acrossFolder = validFolder(fullfile(baseFolder, 'Data', 'AcrossOffsetSummaries'));

pooledPath = fullfile(acrossFolder, 'BetaNoiseRegressionFit.mat');
sessionFitPath = fullfile(acrossFolder, 'BetaNoiseSessionFits.mat');
analysisPath = fullfile(acrossFolder, 'BetaSessionFitAnalysis.mat');

P = load(pooledPath, 'noiseFit');
F = load(sessionFitPath, 'sessionFitTable');
A = load(analysisPath, 'sessionBetaAnalysis');

noiseFit = P.noiseFit;
observedTable = F.sessionFitTable;
observedAnalysis = A.sessionBetaAnalysis;

betaCommon = double(noiseFit.betaNoise);
pooledNames = string(noiseFit.sessionFileNames(:));
pooledAlpha = double(noiseFit.alpha(:));

% Use exactly the sessions included in the session-fit analysis.
validObserved = isfinite(observedTable.beta) & ...
  isfinite(observedTable.betaSE) & ...
  observedTable.betaSE > 0 & ...
  observedTable.exitflag > 0;

observedTable = observedTable(validObserved,:);
sessionNames = string(observedTable.fileName);
nSessions = height(observedTable);

sessionX = cell(nSessions,1);
sessionAlpha = nan(nSessions,1);

for s = 1:nSessions
  name = sessionNames(s);

  pooledMatch = find(strcmpi(pooledNames, name));
  if numel(pooledMatch) ~= 1
    error('simulateCommonBetaSessionHeterogeneity:NameMatch', ...
      '%s matched %d pooled-fit sessions.', name, numel(pooledMatch));
  end

  sessionAlpha(s) = pooledAlpha(pooledMatch);

  D = load(fullfile(sessionFolder, name), 'sessionNoise');
  N = D.sessionNoise;

  use = logical(N.hasPreferredNoise(:));
  x = double(N.trialAnalysis.effectiveNoisePC(use));
  x = x(:);

  if isempty(x) || any(~isfinite(x))
    error('simulateCommonBetaSessionHeterogeneity:BadPredictor', ...
      '%s has invalid preferred-noise predictor values.', name);
  end

  sessionX{s} = x;
end

observedSD = std(observedTable.beta);
observedTau = observedAnalysis.randomEffects.tau;
observedI2 = observedAnalysis.randomEffects.I2;

simMeanBeta = nan(nSim,1);
simMedianBeta = nan(nSim,1);
simSDBeta = nan(nSim,1);
simTau = nan(nSim,1);
simI2 = nan(nSim,1);
simNValid = nan(nSim,1);

for iSim = 1:nSim
  if iSim == 1 || mod(iSim,25) == 0
    fprintf('Simulation %d of %d\n', iSim, nSim);
  end

  betaHat = nan(nSessions,1);
  betaSE = nan(nSessions,1);

  for s = 1:nSessions
    x = sessionX{s};
    eta = sessionAlpha(s) + betaCommon*x;
    p = psychometricProbability(eta);
    y = double(rand(size(p)) < p);

    % Rare all-correct/all-error samples are uninformative.
    if all(y == 1) || all(y == 0)
      continue;
    end

    [betaHat(s), betaSE(s)] = fitOneSession(x, y);
  end

  valid = isfinite(betaHat) & isfinite(betaSE) & betaSE > 0;
  simNValid(iSim) = sum(valid);

  if sum(valid) < 3
    continue;
  end

  b = betaHat(valid);
  se = betaSE(valid);

  simMeanBeta(iSim) = mean(b);
  simMedianBeta(iSim) = median(b);
  simSDBeta(iSim) = std(b);

  [simTau(iSim), simI2(iSim)] = heterogeneityStats(b, se);
end

commonBetaSim = struct();
commonBetaSim.version = 1;
commonBetaSim.nSim = nSim;
commonBetaSim.seed = seed;
commonBetaSim.betaCommon = betaCommon;
commonBetaSim.nullModel = ...
  ['common pooled noise beta with observed effective-noise predictors ' ...
   'and pooled-fit session intercepts'];

commonBetaSim.sessionFileNames = cellstr(sessionNames);
commonBetaSim.sessionAlpha = sessionAlpha;

commonBetaSim.observedSD = observedSD;
commonBetaSim.observedTau = observedTau;
commonBetaSim.observedI2 = observedI2;

commonBetaSim.simMeanBeta = simMeanBeta;
commonBetaSim.simMedianBeta = simMedianBeta;
commonBetaSim.simSDBeta = simSDBeta;
commonBetaSim.simTau = simTau;
commonBetaSim.simI2 = simI2;
commonBetaSim.simNValid = simNValid;

commonBetaSim.sdNullCI95 = prctile(simSDBeta, [2.5 97.5]);
commonBetaSim.tauNullCI95 = prctile(simTau, [2.5 97.5]);
commonBetaSim.i2NullCI95 = prctile(simI2, [2.5 97.5]);

commonBetaSim.pSD = mean(simSDBeta >= observedSD, 'omitnan');
commonBetaSim.pTau = mean(simTau >= observedTau, 'omitnan');
commonBetaSim.pI2 = mean(simI2 >= observedI2, 'omitnan');

commonBetaSim.createdBy = mfilename;
commonBetaSim.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, 'BetaCommonSlopeSimulation.mat');
save(outputPath, 'commonBetaSim', '-v7.3');

fprintf('\nCommon beta used: %.6g\n', betaCommon);
fprintf('Observed session-beta SD: %.6g\n', observedSD);
fprintf('Null SD 95%% interval: %.6g to %.6g; p = %.4g\n', ...
  commonBetaSim.sdNullCI95(1), commonBetaSim.sdNullCI95(2), ...
  commonBetaSim.pSD);

fprintf('Observed tau: %.6g\n', observedTau);
fprintf('Null tau 95%% interval: %.6g to %.6g; p = %.4g\n', ...
  commonBetaSim.tauNullCI95(1), commonBetaSim.tauNullCI95(2), ...
  commonBetaSim.pTau);

fprintf('Observed I^2: %.6g\n', observedI2);
fprintf('Null I^2 95%% interval: %.6g to %.6g; p = %.4g\n', ...
  commonBetaSim.i2NullCI95(1), commonBetaSim.i2NullCI95(2), ...
  commonBetaSim.pI2);

fprintf('Saved %s\n', outputPath);

plotNullDistribution(simSDBeta, observedSD, ...
  'Session beta SD', 'Null distribution of session-beta SD');

plotNullDistribution(simTau, observedTau, ...
  'Random-effects \tau', 'Null distribution of between-session SD');

plotNullDistribution(simI2, observedI2, ...
  'I^2', 'Null distribution of heterogeneity I^2');
end

% -------------------------------------------------------------------------
function [beta, betaSE] = fitOneSession(x, y)

xCenter = mean(x);
xC = x - xCenter;

pObs = mean(y);
q = min(max(2*pObs - 1, 1e-4), 1 - 1e-4);
theta0 = [0.1; log(q/(1-q))];

objective = @(theta) nllGrad(theta, xC, y);

options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', 'off', ...
  'MaxFunctionEvaluations', 2e4, ...
  'MaxIterations', 500, ...
  'FunctionTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-7, ...
  'StepTolerance', 1e-10);

[theta, ~, exitflag, ~, ~, hessian] = ...
  fminunc(objective, theta0, options);

if exitflag <= 0 || any(~isfinite(theta))
  beta = NaN;
  betaSE = NaN;
  return;
end

hessian = full(double(hessian));
covTheta = full(inv(hessian));

beta = full(double(theta(1)));
betaSE = sqrt(covTheta(1,1));

if ~isfinite(betaSE) || betaSE <= 0
  beta = NaN;
  betaSE = NaN;
end
end

% -------------------------------------------------------------------------
function [tau, I2] = heterogeneityStats(beta, betaSE)

v = betaSE.^2;
w = 1 ./ v;

fixedMean = sum(w.*beta) / sum(w);
Q = sum(w.*(beta-fixedMean).^2);
df = numel(beta)-1;

C = sum(w) - sum(w.^2)/sum(w);
tau2 = max(0, (Q-df)/C);

tau = sqrt(tau2);

if Q > 0
  I2 = max(0, (Q-df)/Q);
else
  I2 = 0;
end
end

% -------------------------------------------------------------------------
function [nll, gradient] = nllGrad(theta, x, y)

beta = theta(1);
alpha = theta(2);

eta = alpha + beta*x;
[s, p] = logisticAndPsychometricProbability(eta);

pSafe = min(max(p, realmin), 1-eps);
nll = -sum(y.*log(pSafe) + (1-y).*log(1-pSafe));

dp = 0.5*s.*(1-s);
dEta = (p-y).*dp ./ (pSafe.*(1-pSafe));

gradient = [sum(dEta.*x); sum(dEta)];
end

% -------------------------------------------------------------------------
function p = psychometricProbability(eta)
[~, p] = logisticAndPsychometricProbability(eta);
end

% -------------------------------------------------------------------------
function [s, p] = logisticAndPsychometricProbability(eta)

s = zeros(size(eta));

positive = eta >= 0;
s(positive) = 1 ./ (1 + exp(-eta(positive)));

z = exp(eta(~positive));
s(~positive) = z ./ (1 + z);

p = 0.5 + 0.5*s;
end

% -------------------------------------------------------------------------
function plotNullDistribution(values, observed, xLabelText, titleText)

values = values(isfinite(values));

figure;
histogram(values);
hold on
xline(observed, 'r--', 'Observed', 'LineWidth', 1.5);
xlabel(xLabelText);
ylabel('Number of simulations');
title(titleText);
box off
end
