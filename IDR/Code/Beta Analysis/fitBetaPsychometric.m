function psychometricFit = fitBetaPsychometric()
% fitBetaPsychometric
% Fit a pooled 2AFC psychometric function with analytic gradient.
%
% Model:
%   P(correct) = 0.5 + 0.5 ./ (1 + exp(-(alpha_session + beta*x)))
%
% x is trialwise effective coherence. The fit has one shared slope beta and
% one intercept per session. The 0.5 chance floor is fixed.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%
% Optionally uses an existing BetaPsychometricFit.mat as the starting point.
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaPsychometricFit.mat

cleanupObj = initProjectPath(); %#ok<NASGU>

baseFolder = folderPath();
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
outputFolder = validFolder(fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions'));
outputPath = fullfile(outputFolder, 'BetaPsychometricFit.mat');

files = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(files)
  error('fitBetaPsychometric:NoSessionFiles', ...
    'No SessionData files found in %s.', sessionDataFolder);
end

[~, order] = sort({files.name});
files = files(order);
nSessions = numel(files);

effectiveCohPC = [];
correct = [];
sessionIndex = [];

sessionNTrials = zeros(nSessions, 1);
sessionNCorrect = zeros(nSessions, 1);
sessionMeanSignalCohPC = nan(nSessions, 1);
sessionSignalCoherencesPC = cell(nSessions, 1);

for iSession = 1:nSessions
  S = load(fullfile(sessionDataFolder, files(iSession).name), 'sessionNoise');
  if ~isfield(S, 'sessionNoise')
    error('fitBetaPsychometric:MissingSessionNoise', ...
      '%s does not contain sessionNoise.', files(iSession).name);
  end

  N = S.sessionNoise;
  if ~isfield(N, 'trialAnalysis') || ...
      ~isfield(N.trialAnalysis, 'effectiveCohPC')
    error('fitBetaPsychometric:MissingEffectiveCoherence', ...
      '%s lacks trialAnalysis.effectiveCohPC.', files(iSession).name);
  end

  x = double(N.trialAnalysis.effectiveCohPC(:));
  y = double(N.trialOutcome(:) == 0);

  if numel(x) ~= numel(y) || any(~isfinite(x))
    error('fitBetaPsychometric:BadTrialData', ...
      '%s contains invalid trial data.', files(iSession).name);
  end

  n = numel(y);
  effectiveCohPC = [effectiveCohPC; x]; %#ok<AGROW>
  correct = [correct; y]; %#ok<AGROW>
  sessionIndex = [sessionIndex; repmat(iSession, n, 1)]; %#ok<AGROW>

  sessionNTrials(iSession) = n;
  sessionNCorrect(iSession) = sum(y);
  sessionMeanSignalCohPC(iSession) = mean(double(N.signalCohPC(:)));
  sessionSignalCoherencesPC{iSession} = unique(double(N.signalCohPC(:)))';
end

nTrials = numel(correct);

% Center predictor for numerical stability.
xCenter = mean(effectiveCohPC);
xCentered = effectiveCohPC - xCenter;

% Prefer the previous fit as the starting point.
theta0 = [];
if isfile(outputPath)
  old = load(outputPath, 'psychometricFit');
  if isfield(old, 'psychometricFit')
    F0 = old.psychometricFit;
    if isfield(F0, 'beta') && isfield(F0, 'alpha') && ...
        numel(F0.alpha) == nSessions
      beta0 = double(F0.beta);
      alphaCentered0 = double(F0.alpha(:)) + beta0 * xCenter;
      theta0 = [beta0; alphaCentered0];
      fprintf('Starting from the previous saved fit.\n');
    end
  end
end

if isempty(theta0)
  beta0 = 0.1;
  alphaCentered0 = nan(nSessions, 1);

  for iSession = 1:nSessions
    p = mean(correct(sessionIndex == iSession));
    q = min(max(2*p - 1, 1e-4), 1 - 1e-4);
    alphaCentered0(iSession) = log(q / (1 - q));
  end

  theta0 = [beta0; alphaCentered0];
end

objective = @(theta) negativeLogLikelihoodAndGradient( ...
  theta, xCentered, correct, sessionIndex, nSessions);

options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', 'iter', ...
  'MaxFunctionEvaluations', 5e5, ...
  'MaxIterations', 1e4, ...
  'OptimalityTolerance', 1e-9, ...
  'StepTolerance', 1e-12);

[thetaHat, nll, exitflag, output, gradient, hessian] = ...
  fminunc(objective, theta0, options);

gradient = full(gradient);
hessian = full(hessian);


beta = thetaHat(1);
alphaCentered = thetaHat(2:end);
alpha = alphaCentered - beta*xCenter;
threshold75PC = -alpha / beta;

% Covariance from observed Hessian.
covTheta = inv(hessian);
seTheta = sqrt(diag(covTheta));
betaSE = seTheta(1);

J = eye(nSessions + 1);
J(2:end, 1) = -xCenter;
covOriginal = J * covTheta * J';

alphaSE = sqrt(diag(covOriginal(2:end, 2:end)));

thresholdSE = nan(nSessions, 1);
for iSession = 1:nSessions
  a = alpha(iSession);
  g = zeros(nSessions + 1, 1);
  g(1) = a / beta^2;
  g(iSession + 1) = -1 / beta;
  thresholdSE(iSession) = sqrt(g' * covOriginal * g);
end

betaCI95 = beta + 1.96*betaSE*[-1 1];
alphaCI95 = [alpha - 1.96*alphaSE, alpha + 1.96*alphaSE];
thresholdCI95 = [threshold75PC - 1.96*thresholdSE, ...
                 threshold75PC + 1.96*thresholdSE];

eta = alpha(sessionIndex) + beta*effectiveCohPC;
fittedProbability = psychometricProbability(eta);

logLikelihood = -nll;
nParameters = nSessions + 1;
AIC = 2*nParameters - 2*logLikelihood;
BIC = log(nTrials)*nParameters - 2*logLikelihood;

sessionObservedPCorrect = sessionNCorrect ./ sessionNTrials;
sessionPredictedPCorrect = accumarray(sessionIndex, fittedProbability, ...
  [nSessions 1], @mean);

beta = full(double(beta));
betaSE = full(double(betaSE));
betaCI95 = full(double(betaCI95));

alpha = full(double(alpha));
alphaSE = full(double(alphaSE));
alphaCI95 = full(double(alphaCI95));

threshold75PC = full(double(threshold75PC));
thresholdSE = full(double(thresholdSE));
thresholdCI95 = full(double(thresholdCI95));

gradient = full(double(gradient));
hessian = full(double(hessian));
covOriginal = full(double(covOriginal));

psychometricFit = struct();
psychometricFit.version = 2;
psychometricFit.model = ...
  'P(correct)=0.5+0.5*logistic(alpha_session+beta*effectiveCohPC)';
psychometricFit.chanceFloor = 0.5;

psychometricFit.beta = beta;
psychometricFit.betaSE = betaSE;
psychometricFit.betaCI95 = betaCI95;
psychometricFit.alpha = alpha;
psychometricFit.alphaSE = alphaSE;
psychometricFit.alphaCI95 = alphaCI95;
psychometricFit.threshold75PC = threshold75PC;
psychometricFit.threshold75SE = thresholdSE;
psychometricFit.threshold75CI95 = thresholdCI95;

psychometricFit.sessionFileNames = {files.name};
psychometricFit.sessionNTrials = sessionNTrials;
psychometricFit.sessionNCorrect = sessionNCorrect;
psychometricFit.sessionObservedPCorrect = sessionObservedPCorrect;
psychometricFit.sessionPredictedPCorrect = sessionPredictedPCorrect;
psychometricFit.sessionMeanSignalCohPC = sessionMeanSignalCohPC;
psychometricFit.sessionSignalCoherencesPC = sessionSignalCoherencesPC;

psychometricFit.effectiveCohPC = effectiveCohPC;
psychometricFit.correct = correct;
psychometricFit.sessionIndex = sessionIndex;
psychometricFit.fittedProbability = fittedProbability;

psychometricFit.nTrials = nTrials;
psychometricFit.nSessions = nSessions;
psychometricFit.logLikelihood = logLikelihood;
psychometricFit.AIC = AIC;
psychometricFit.BIC = BIC;
psychometricFit.exitflag = exitflag;
psychometricFit.optimizerOutput = output;
psychometricFit.gradient = gradient;
psychometricFit.gradientInfNorm = norm(gradient, Inf);
psychometricFit.hessian = hessian;
psychometricFit.covariance = covOriginal;
psychometricFit.xCenter = xCenter;
psychometricFit.createdBy = mfilename;
psychometricFit.createdDate = datetime('now');

fprintf('Saved %s\n', outputPath);
fprintf('Shared beta: %.6g (SE %.6g; 95%% CI %.6g to %.6g)\n', ...
  beta, betaSE, betaCI95(1), betaCI95(2));
fprintf('Gradient infinity norm: %.6g\n', psychometricFit.gradientInfNorm);
fprintf('Exitflag: %d\n', exitflag);

save(outputPath, 'psychometricFit', '-v7.3');

plotPsychometricFit(psychometricFit);
end

% -------------------------------------------------------------------------
function [nll, gradient] = negativeLogLikelihoodAndGradient( ...
  theta, x, y, sessionIndex, nSessions)

beta = theta(1);
alpha = theta(2:end);

eta = alpha(sessionIndex) + beta*x;
[s, p] = logisticAndPsychometricProbability(eta);

pSafe = min(max(p, realmin), 1 - eps);
nll = -sum(y .* log(pSafe) + (1-y) .* log(1-pSafe));

% dp/deta for p = 0.5 + 0.5*s
dp = 0.5 * s .* (1 - s);

% d(NLL)/deta
dEta = (p - y) .* dp ./ (pSafe .* (1 - pSafe));

gradient = zeros(nSessions + 1, 1);
gradient(1) = sum(dEta .* x);
gradient(2:end) = accumarray(sessionIndex, dEta, ...
  [nSessions 1], @sum, 0);
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
function plotPsychometricFit(F)

thresholdByTrial = F.threshold75PC(F.sessionIndex);
centeredCohPC = F.effectiveCohPC - thresholdByTrial;

nBins = 15;
edges = linspace(prctile(centeredCohPC, 1), ...
                 prctile(centeredCohPC, 99), nBins + 1);

bin = discretize(centeredCohPC, edges);
binCenter = nan(nBins, 1);
binPCorrect = nan(nBins, 1);
binN = zeros(nBins, 1);

for iBin = 1:nBins
  use = bin == iBin;
  if any(use)
    binCenter(iBin) = mean(centeredCohPC(use));
    binPCorrect(iBin) = mean(F.correct(use));
    binN(iBin) = sum(use);
  end
end

xGrid = linspace(min(edges), max(edges), 400);
pGrid = psychometricProbability(F.beta*xGrid);

figure;
hold on
plot(xGrid, pGrid, 'k-', 'LineWidth', 1.5);
plot(binCenter, binPCorrect, 'ko', ...
  'MarkerFaceColor', 'w', 'LineWidth', 1);
yline(0.5, 'k:');
yline(0.75, 'k--');
xline(0, 'k:');

xlabel('Effective coherence relative to session threshold (%)');
ylabel('Proportion correct');
title(sprintf('Shared psychometric slope = %.4g', F.beta));
ylim([0.45 1.02]);
box off

for iBin = 1:nBins
  if isfinite(binCenter(iBin))
    text(binCenter(iBin), binPCorrect(iBin)-0.025, ...
      sprintf('%d', binN(iBin)), ...
      'HorizontalAlignment', 'center', ...
      'VerticalAlignment', 'top', ...
      'FontSize', 7);
  end
end
end
