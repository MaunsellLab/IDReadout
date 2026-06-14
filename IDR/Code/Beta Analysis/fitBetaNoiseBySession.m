function sessionFits = fitBetaNoiseBySession()
% fitBetaNoiseBySession
% Fit one preferred-noise regression beta per session.
%
% Each session fit uses:
%   - only trials with hasPreferredNoise == true
%   - that session's stored leave-one-out effectiveNoisePC
%   - the same fixed 0.5 floor used in the pooled fits
%
% Model for one session:
%   P(correct) = 0.5 + 0.5*logistic(alpha + beta*effectiveNoisePC)
%
% Each session beta is compared with the single pooled psychometric beta
% from BetaPsychometricFit.mat. That reference beta is not recomputed with
% the current session omitted.
%
% Reads:
%   Data/FullSessions/BetaAnalysis/SessionData/*.mat
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaPsychometricFit.mat
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaNoiseSessionFits.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

baseFolder = domainFolder(mfilename('fullpath'));

sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');

acrossFolder = validFolder(fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions'));

psychometricPath = fullfile(acrossFolder, 'BetaPsychometricFit.mat');
if ~isfile(psychometricPath)
  error('fitBetaNoiseBySession:MissingReferenceFit', ...
    'Psychometric fit not found: %s', psychometricPath);
end

P = load(psychometricPath, 'psychometricFit');
if ~isfield(P, 'psychometricFit')
  error('fitBetaNoiseBySession:MissingPsychometricFit', ...
    '%s does not contain psychometricFit.', psychometricPath);
end

referenceBeta = double(P.psychometricFit.beta);

files = dir(fullfile(sessionDataFolder, '*.mat'));
if isempty(files)
  error('fitBetaNoiseBySession:NoSessionFiles', ...
    'No SessionData files found in %s.', sessionDataFolder);
end

[~, order] = sort({files.name});
files = files(order);
nSessions = numel(files);

sessionFits = repmat(emptyFitRow(), nSessions, 1);

for iSession = 1:nSessions
  fileName = files(iSession).name;
  filePath = fullfile(sessionDataFolder, fileName);

  S = load(filePath, 'sessionNoise');
  if ~isfield(S, 'sessionNoise')
    error('fitBetaNoiseBySession:MissingSessionNoise', ...
      '%s does not contain sessionNoise.', fileName);
  end

  N = S.sessionNoise;

  if ~isfield(N, 'trialAnalysis') || ...
      ~isfield(N.trialAnalysis, 'effectiveNoisePC')
    error('fitBetaNoiseBySession:MissingEffectiveNoise', ...
      '%s lacks trialAnalysis.effectiveNoisePC.', fileName);
  end

  useTrials = logical(N.hasPreferredNoise(:));

  x = double(N.trialAnalysis.effectiveNoisePC(useTrials));
  y = double(N.trialOutcome(useTrials) == 0);

  x = x(:);
  y = y(:);

  if isempty(x)
    error('fitBetaNoiseBySession:NoPreferredNoiseTrials', ...
      '%s contains no preferred-noise trials.', fileName);
  end

  if numel(x) ~= numel(y) || any(~isfinite(x)) || ...
      any(~ismember(y, [0 1]))
    error('fitBetaNoiseBySession:BadTrialData', ...
      '%s contains invalid session-fit data.', fileName);
  end

  nTrials = numel(y);
  nCorrect = sum(y);
  nError = nTrials - nCorrect;

  if nCorrect == 0 || nError == 0
    warning('fitBetaNoiseBySession:EmptyOutcomeClass', ...
      '%s has %d correct and %d error trials; returning NaNs.', ...
      fileName, nCorrect, nError);

    row = emptyFitRow();
    row.fileName = fileName;
    row.nTrials = nTrials;
    row.nCorrect = nCorrect;
    row.nError = nError;
    row.fractionCorrect = mean(y);
    row.referenceBeta = referenceBeta;
    sessionFits(iSession) = row;
    continue;
  end

  [beta, alpha, covTheta, exitflag, output, gradient, hessian] = ...
    fitOneSession(x, y);

  betaSE = sqrt(covTheta(1,1));
  alphaSE = sqrt(covTheta(2,2));

  betaCI95 = beta + 1.96*betaSE*[-1 1];
  alphaCI95 = alpha + 1.96*alphaSE*[-1 1];

  eta = alpha + beta*x;
  fittedProbability = psychometricProbability(eta);

  threshold75PC = -alpha / beta;

  g = [alpha / beta^2; -1 / beta];
  thresholdSE = sqrt(g' * covTheta * g);
  thresholdCI95 = threshold75PC + 1.96*thresholdSE*[-1 1];

  row = emptyFitRow();
  row.fileName = fileName;
  row.nTrials = nTrials;
  row.nCorrect = nCorrect;
  row.nError = nError;
  row.fractionCorrect = mean(y);

  row.beta = beta;
  row.betaSE = betaSE;
  row.betaCI95Low = betaCI95(1);
  row.betaCI95High = betaCI95(2);

  row.alpha = alpha;
  row.alphaSE = alphaSE;
  row.alphaCI95Low = alphaCI95(1);
  row.alphaCI95High = alphaCI95(2);

  row.threshold75PC = threshold75PC;
  row.threshold75SE = thresholdSE;
  row.threshold75CI95Low = thresholdCI95(1);
  row.threshold75CI95High = thresholdCI95(2);

  row.referenceBeta = referenceBeta;
  row.betaDifference = beta - referenceBeta;
  row.betaRatio = beta / referenceBeta;

  row.meanEffectiveNoisePC = mean(x);
  row.sdEffectiveNoisePC = std(x);
  row.meanFittedProbability = mean(fittedProbability);

  row.exitflag = exitflag;
  row.firstOrderOptimality = output.firstorderopt;
  row.gradientInfNorm = norm(gradient, Inf);
  row.iterations = output.iterations;

  sessionFits(iSession) = row;
end

sessionFitTable = struct2table(sessionFits);

summary = struct();
summary.version = 1;
summary.model = ...
  'P(correct)=0.5+0.5*logistic(alpha+beta*effectiveNoisePC)';
summary.referenceBeta = referenceBeta;
summary.referenceSource = psychometricPath;
summary.sessionFitTable = sessionFitTable;

validBeta = isfinite(sessionFitTable.beta);
summary.nSessions = height(sessionFitTable);
summary.nValidSessions = sum(validBeta);
summary.meanBeta = mean(sessionFitTable.beta(validBeta));
summary.sdBeta = std(sessionFitTable.beta(validBeta));
summary.semBeta = summary.sdBeta / sqrt(summary.nValidSessions);
summary.meanDifference = mean(sessionFitTable.betaDifference(validBeta));
summary.meanRatio = mean(sessionFitTable.betaRatio(validBeta));
summary.createdBy = mfilename;
summary.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, 'BetaNoiseSessionFits.mat');
save(outputPath, 'summary', 'sessionFitTable', '-v7.3');

fprintf('Saved %s\n', outputPath);
fprintf('Reference psychometric beta: %.6g\n', referenceBeta);
fprintf('Valid session fits: %d of %d\n', ...
  summary.nValidSessions, summary.nSessions);
fprintf('Mean session beta: %.6g (SD %.6g; SEM %.6g)\n', ...
  summary.meanBeta, summary.sdBeta, summary.semBeta);
fprintf('Mean beta/reference ratio: %.6g\n', summary.meanRatio);

plotSessionFits(sessionFitTable, referenceBeta);
end

% -------------------------------------------------------------------------
function [beta, alpha, covTheta, exitflag, output, gradient, hessian] = ...
  fitOneSession(x, y)

xCenter = mean(x);
xCentered = x - xCenter;

p = mean(y);
q = min(max(2*p - 1, 1e-4), 1 - 1e-4);

theta0 = [0.1; log(q/(1-q))];

objective = @(theta) nllGrad(theta, xCentered, y);

options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', 'off', ...
  'MaxFunctionEvaluations', 5e4, ...
  'MaxIterations', 1000, ...
  'FunctionTolerance', 1e-12, ...
  'OptimalityTolerance', 1e-8, ...
  'StepTolerance', 1e-12);

[theta, ~, exitflag, output, gradient, hessian] = ...
  fminunc(objective, theta0, options);

gradient = full(double(gradient));
hessian = full(double(hessian));
covCentered = full(inv(hessian));

beta = full(double(theta(1)));
alphaCentered = full(double(theta(2)));
alpha = alphaCentered - beta*xCenter;

J = [1 0; -xCenter 1];
covTheta = J * covCentered * J';
covTheta = full(double(covTheta));
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
function row = emptyFitRow()

row = struct( ...
  'fileName', '', ...
  'nTrials', NaN, ...
  'nCorrect', NaN, ...
  'nError', NaN, ...
  'fractionCorrect', NaN, ...
  'beta', NaN, ...
  'betaSE', NaN, ...
  'betaCI95Low', NaN, ...
  'betaCI95High', NaN, ...
  'alpha', NaN, ...
  'alphaSE', NaN, ...
  'alphaCI95Low', NaN, ...
  'alphaCI95High', NaN, ...
  'threshold75PC', NaN, ...
  'threshold75SE', NaN, ...
  'threshold75CI95Low', NaN, ...
  'threshold75CI95High', NaN, ...
  'referenceBeta', NaN, ...
  'betaDifference', NaN, ...
  'betaRatio', NaN, ...
  'meanEffectiveNoisePC', NaN, ...
  'sdEffectiveNoisePC', NaN, ...
  'meanFittedProbability', NaN, ...
  'exitflag', NaN, ...
  'firstOrderOptimality', NaN, ...
  'gradientInfNorm', NaN, ...
  'iterations', NaN);
end

% -------------------------------------------------------------------------
function plotSessionFits(T, referenceBeta)

valid = isfinite(T.beta) & isfinite(T.betaSE);
x = find(valid);

figure;
hold on

for i = 1:numel(x)
  k = x(i);
  plot([k k], ...
    [T.betaCI95Low(k), T.betaCI95High(k)], ...
    'k-');
end

plot(x, T.beta(valid), 'ko', ...
  'MarkerFaceColor', 'w', ...
  'MarkerSize', 4);

yline(referenceBeta, 'r--', ...
  sprintf('Reference \\beta = %.4g', referenceBeta));

xlabel('Session');
ylabel('Noise regression beta');
title(sprintf('Session noise betas (%d sessions)', sum(valid)));

box off
end
