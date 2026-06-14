function noiseFit = fitBetaNoiseRegression()
% fitBetaNoiseRegression
% Fit the trialwise preferred-noise regression beta.
%
% Model:
%   P(correct) = 0.5 + 0.5 ./ ...
%       (1 + exp(-(alpha_session + betaNoise*effectiveNoisePC)))
%
% Uses only trials with:
%   sessionNoise.hasPreferredNoise == true
%
% The model has one shared noise slope and one intercept per session.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaNoiseRegressionFit.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

baseFolder = domainFolder(mfilename('fullpath'));
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
outputFolder = validFolder(fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions'));
outputPath = fullfile(outputFolder, 'BetaNoiseRegressionFit.mat');

files = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(files)
  error('fitBetaNoiseRegression:NoSessionFiles', ...
    'No SessionData files found in %s.', sessionDataFolder);
end

[~, order] = sort({files.name});
files = files(order);
nSessions = numel(files);

effectiveNoisePC = [];
correct = [];
sessionIndex = [];

sessionNTrials = zeros(nSessions, 1);
sessionNCorrect = zeros(nSessions, 1);
sessionSignalCohPC = cell(nSessions, 1);

for iSession = 1:nSessions
  filePath = fullfile(sessionDataFolder, files(iSession).name);
  S = load(filePath, 'sessionNoise');

  if ~isfield(S, 'sessionNoise')
    error('fitBetaNoiseRegression:MissingSessionNoise', ...
      '%s does not contain sessionNoise.', files(iSession).name);
  end

  N = S.sessionNoise;

  if ~isfield(N, 'trialAnalysis') || ...
      ~isfield(N.trialAnalysis, 'effectiveNoisePC')
    error('fitBetaNoiseRegression:MissingEffectiveNoise', ...
      '%s lacks trialAnalysis.effectiveNoisePC.', files(iSession).name);
  end

  useTrials = logical(N.hasPreferredNoise(:));

  x = double(N.trialAnalysis.effectiveNoisePC(useTrials));
  y = double(N.trialOutcome(useTrials) == 0);
  signal = double(N.signalCohPC(useTrials));

  x = x(:);
  y = y(:);
  signal = signal(:);

  if isempty(x)
    error('fitBetaNoiseRegression:NoPreferredNoiseTrials', ...
      '%s contains no preferred-noise trials.', files(iSession).name);
  end

  if numel(x) ~= numel(y) || any(~isfinite(x)) || ...
      any(~ismember(y, [0 1]))
    error('fitBetaNoiseRegression:BadTrialData', ...
      '%s contains invalid noise-regression data.', files(iSession).name);
  end

  % Signal coherence should be constant among the preferred-noise trials
  % within each session. The session intercept then absorbs that operating
  % point.
  signalValues = unique(signal);
  if numel(signalValues) ~= 1
    error('fitBetaNoiseRegression:MultipleNoiseSignalCoherences', ...
      ['%s has %d signal coherences among preferred-noise trials; ' ...
       'expected exactly one.'], ...
      files(iSession).name, numel(signalValues));
  end

  n = numel(y);

  effectiveNoisePC = [effectiveNoisePC; x]; %#ok<AGROW>
  correct = [correct; y]; %#ok<AGROW>
  sessionIndex = [sessionIndex; repmat(iSession, n, 1)]; %#ok<AGROW>

  sessionNTrials(iSession) = n;
  sessionNCorrect(iSession) = sum(y);
  sessionSignalCohPC{iSession} = signalValues(:)';
end

nTrials = numel(correct);

% Center predictor for numerical stability.
xCenter = mean(effectiveNoisePC);
xCentered = effectiveNoisePC - xCenter;

theta0 = [];

if isfile(outputPath)
  old = load(outputPath, 'noiseFit');
  if isfield(old, 'noiseFit')
    F0 = old.noiseFit;
    if isfield(F0, 'betaNoise') && isfield(F0, 'alpha') && ...
        numel(F0.alpha) == nSessions
      beta0 = double(F0.betaNoise);
      alphaCentered0 = double(F0.alpha(:)) + beta0*xCenter;
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
  'FunctionTolerance', 1e-12, ...
  'StepTolerance', 1e-12);

[thetaHat, nll, exitflag, output, gradient, hessian] = ...
  fminunc(objective, theta0, options);

gradient = full(double(gradient));
hessian = full(double(hessian));

betaNoise = full(double(thetaHat(1)));
alphaCentered = full(double(thetaHat(2:end)));
alpha = alphaCentered - betaNoise*xCenter;

covTheta = inv(hessian);
covTheta = full(double(covTheta));

seTheta = sqrt(diag(covTheta));
betaNoiseSE = full(double(seTheta(1)));

J = eye(nSessions + 1);
J(2:end, 1) = -xCenter;
covOriginal = J * covTheta * J';
covOriginal = full(double(covOriginal));

alphaSE = sqrt(diag(covOriginal(2:end, 2:end)));

betaNoiseCI95 = betaNoise + 1.96*betaNoiseSE*[-1 1];
alphaCI95 = [alpha - 1.96*alphaSE, alpha + 1.96*alphaSE];

eta = alpha(sessionIndex) + betaNoise*effectiveNoisePC;
fittedProbability = psychometricProbability(eta);

logLikelihood = -nll;
nParameters = nSessions + 1;
AIC = 2*nParameters - 2*logLikelihood;
BIC = log(nTrials)*nParameters - 2*logLikelihood;

sessionObservedPCorrect = sessionNCorrect ./ sessionNTrials;
sessionPredictedPCorrect = accumarray(sessionIndex, fittedProbability, ...
  [nSessions 1], @mean);

noiseFit = struct();
noiseFit.version = 1;
noiseFit.model = ...
  'P(correct)=0.5+0.5*logistic(alpha_session+betaNoise*effectiveNoisePC)';
noiseFit.chanceFloor = 0.5;

noiseFit.betaNoise = betaNoise;
noiseFit.betaNoiseSE = betaNoiseSE;
noiseFit.betaNoiseCI95 = betaNoiseCI95;

noiseFit.alpha = alpha;
noiseFit.alphaSE = alphaSE;
noiseFit.alphaCI95 = alphaCI95;

noiseFit.sessionFileNames = {files.name};
noiseFit.sessionNTrials = sessionNTrials;
noiseFit.sessionNCorrect = sessionNCorrect;
noiseFit.sessionObservedPCorrect = sessionObservedPCorrect;
noiseFit.sessionPredictedPCorrect = sessionPredictedPCorrect;
noiseFit.sessionSignalCohPC = sessionSignalCohPC;

noiseFit.effectiveNoisePC = effectiveNoisePC;
noiseFit.correct = correct;
noiseFit.sessionIndex = sessionIndex;
noiseFit.fittedProbability = fittedProbability;

noiseFit.nTrials = nTrials;
noiseFit.nSessions = nSessions;
noiseFit.logLikelihood = logLikelihood;
noiseFit.AIC = AIC;
noiseFit.BIC = BIC;

noiseFit.exitflag = exitflag;
noiseFit.optimizerOutput = output;
noiseFit.gradient = gradient;
noiseFit.gradientInfNorm = norm(gradient, Inf);
noiseFit.hessian = hessian;
noiseFit.covariance = covOriginal;
noiseFit.xCenter = xCenter;

noiseFit.createdBy = mfilename;
noiseFit.createdDate = datetime('now');

fprintf('Shared noise beta: %.6g (SE %.6g; 95%% CI %.6g to %.6g)\n', ...
  betaNoise, betaNoiseSE, betaNoiseCI95(1), betaNoiseCI95(2));
fprintf('Gradient infinity norm: %.6g\n', noiseFit.gradientInfNorm);
fprintf('Exitflag: %d\n', exitflag);
fprintf('Included %d preferred-noise trials from %d sessions.\n', ...
  nTrials, nSessions);

save(outputPath, 'noiseFit', '-v7.3');
fprintf('Saved %s\n', outputPath);

plotNoiseRegression(noiseFit);
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

dp = 0.5*s.*(1-s);
dEta = (p-y).*dp ./ (pSafe.*(1-pSafe));

gradient = zeros(nSessions + 1, 1);
gradient(1) = sum(dEta.*x);
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
function plotNoiseRegression(F)

nBins = 15;
edges = linspace(prctile(F.effectiveNoisePC, 1), ...
                 prctile(F.effectiveNoisePC, 99), nBins + 1);

bin = discretize(F.effectiveNoisePC, edges);
binCenter = nan(nBins, 1);
binPCorrect = nan(nBins, 1);
binN = zeros(nBins, 1);

for iBin = 1:nBins
  use = bin == iBin;

  if any(use)
    binCenter(iBin) = mean(F.effectiveNoisePC(use));
    binPCorrect(iBin) = mean(F.correct(use));
    binN(iBin) = sum(use);
  end
end

% Display a curve using the mean fitted session intercept. This is only a
% visualization; the fitted model retains all session-specific intercepts.
meanAlpha = mean(F.alpha);
xGrid = linspace(min(edges), max(edges), 400);
pGrid = psychometricProbability(meanAlpha + F.betaNoise*xGrid);

figure;
hold on

plot(xGrid, pGrid, 'k-', 'LineWidth', 1.5);
plot(binCenter, binPCorrect, 'ko', ...
  'MarkerFaceColor', 'w', 'LineWidth', 1);

xlabel('Effective preferred noise (%)');
ylabel('Proportion correct');
title(sprintf('Shared noise beta = %.4g', F.betaNoise));

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
