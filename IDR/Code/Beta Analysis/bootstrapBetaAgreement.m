function agreement = bootstrapBetaAgreement(nBoot, seed)
% bootstrapBetaAgreement
% Hierarchical paired bootstrap comparing psychometric and noise betas.
%
% Resampling:
%   1. Resample sessions with replacement.
%   2. Within each selected session, resample all increment trials with replacement.
%   3. Fit the psychometric beta to all resampled trials.
%   4. Fit the noise beta to the preferred-noise subset of the same resample.
%
% Both models use:
%   P(correct) = 0.5 + 0.5*logistic(alpha_session + beta*x)
%
% Usage:
%   agreement = bootstrapBetaAgreement();
%   agreement = bootstrapBetaAgreement(1000, 1);
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/BetaAgreementBootstrap.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

if nargin < 1 || isempty(nBoot)
  nBoot = 1000;
end
if nargin < 2 || isempty(seed)
  seed = 1;
end

rng(seed);

baseFolder = domainFolder(mfilename('fullpath'));
sessionDataFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
acrossFolder = validFolder(fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions'));

files = dir(fullfile(sessionDataFolder, '*.mat'));
[~, order] = sort({files.name});
files = files(order);
nSessions = numel(files);

% Load once.
sessionData = cell(nSessions, 1);
for i = 1:nSessions
  S = load(fullfile(sessionDataFolder, files(i).name), 'sessionNoise');
  N = S.sessionNoise;

  D = struct();
  D.correct = double(N.trialOutcome(:) == 0);
  D.effectiveCohPC = double(N.trialAnalysis.effectiveCohPC(:));
  D.effectiveNoisePC = double(N.trialAnalysis.effectiveNoisePC(:));
  D.hasPreferredNoise = logical(N.hasPreferredNoise(:));

  sessionData{i} = D;
end

% Point estimates from saved fits.
P = load(fullfile(acrossFolder, 'BetaPsychometricFit.mat'), 'psychometricFit');
R = load(fullfile(acrossFolder, 'BetaNoiseRegressionFit.mat'), 'noiseFit');

betaPsychPoint = P.psychometricFit.beta;
betaNoisePoint = R.noiseFit.betaNoise;
differencePoint = betaNoisePoint - betaPsychPoint;
ratioPoint = betaNoisePoint / betaPsychPoint;

bootPsych = nan(nBoot,1);
bootNoise = nan(nBoot,1);
bootDifference = nan(nBoot,1);
bootRatio = nan(nBoot,1);
exitPsych = nan(nBoot,1);
exitNoise = nan(nBoot,1);

for b = 1:nBoot
  if b == 1 || mod(b,25) == 0
    fprintf('Bootstrap %d of %d\n', b, nBoot);
  end

  sampledSessions = randi(nSessions, nSessions, 1);

  xPsych = [];
  yPsych = [];
  sPsych = [];

  xNoise = [];
  yNoise = [];
  sNoise = [];

  for j = 1:nSessions
    D = sessionData{sampledSessions(j)};
    n = numel(D.correct);
    idx = randi(n, n, 1);

    y = D.correct(idx);
    xP = D.effectiveCohPC(idx);
    xN = D.effectiveNoisePC(idx);
    hasN = D.hasPreferredNoise(idx);

    xPsych = [xPsych; xP]; %#ok<AGROW>
    yPsych = [yPsych; y]; %#ok<AGROW>
    sPsych = [sPsych; repmat(j,n,1)]; %#ok<AGROW>

    xNoise = [xNoise; xN(hasN)]; %#ok<AGROW>
    yNoise = [yNoise; y(hasN)]; %#ok<AGROW>
    sNoise = [sNoise; repmat(j,sum(hasN),1)]; %#ok<AGROW>
  end

  [bootPsych(b), exitPsych(b)] = fitSharedSlope(xPsych, yPsych, sPsych, nSessions);
  [bootNoise(b), exitNoise(b)] = fitSharedSlope(xNoise, yNoise, sNoise, nSessions);

  bootDifference(b) = bootNoise(b) - bootPsych(b);
  bootRatio(b) = bootNoise(b) / bootPsych(b);
end

agreement = struct();
agreement.version = 1;
agreement.nBoot = nBoot;
agreement.seed = seed;

agreement.betaPsychometric = betaPsychPoint;
agreement.betaNoise = betaNoisePoint;
agreement.difference = differencePoint;
agreement.ratio = ratioPoint;

agreement.bootPsychometric = bootPsych;
agreement.bootNoise = bootNoise;
agreement.bootDifference = bootDifference;
agreement.bootRatio = bootRatio;

agreement.differenceCI95 = prctile(bootDifference, [2.5 97.5]);
agreement.ratioCI95 = prctile(bootRatio, [2.5 97.5]);

agreement.pDifferenceTwoSided = ...
  2 * min(mean(bootDifference <= 0), mean(bootDifference >= 0));

agreement.exitPsychometric = exitPsych;
agreement.exitNoise = exitNoise;
agreement.sessionFileNames = {files.name};
agreement.createdBy = mfilename;
agreement.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, 'BetaAgreementBootstrap.mat');
save(outputPath, 'agreement', '-v7.3');

fprintf('\nPsychometric beta: %.6g\n', betaPsychPoint);
fprintf('Noise beta:        %.6g\n', betaNoisePoint);
fprintf('Difference:        %.6g (95%% CI %.6g to %.6g)\n', ...
  differencePoint, agreement.differenceCI95(1), agreement.differenceCI95(2));
fprintf('Ratio:             %.6g (95%% CI %.6g to %.6g)\n', ...
  ratioPoint, agreement.ratioCI95(1), agreement.ratioCI95(2));
fprintf('Two-sided bootstrap p for difference from zero: %.6g\n', ...
  agreement.pDifferenceTwoSided);
fprintf('Saved %s\n', outputPath);
end

% -------------------------------------------------------------------------
function [beta, exitflag] = fitSharedSlope(x, y, sessionIndex, nSessions)

x = double(x(:));
y = double(y(:));
sessionIndex = double(sessionIndex(:));

xCenter = mean(x);
xC = x - xCenter;

beta0 = 0.1;
alpha0 = nan(nSessions,1);

for s = 1:nSessions
  p = mean(y(sessionIndex == s));
  q = min(max(2*p - 1, 1e-4), 1 - 1e-4);
  alpha0(s) = log(q/(1-q));
end

theta0 = [beta0; alpha0];

objective = @(theta) nllGrad(theta, xC, y, sessionIndex, nSessions);

options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', 'off', ...
  'MaxFunctionEvaluations', 2e5, ...
  'MaxIterations', 2000, ...
  'FunctionTolerance', 1e-10, ...
  'OptimalityTolerance', 1e-7, ...
  'StepTolerance', 1e-10);

[theta, ~, exitflag] = fminunc(objective, theta0, options);
beta = full(double(theta(1)));
end

% -------------------------------------------------------------------------
function [nll, gradient] = nllGrad(theta, x, y, sessionIndex, nSessions)

beta = theta(1);
alpha = theta(2:end);

eta = alpha(sessionIndex) + beta*x;

s = zeros(size(eta));
positive = eta >= 0;
s(positive) = 1 ./ (1 + exp(-eta(positive)));
z = exp(eta(~positive));
s(~positive) = z ./ (1 + z);

p = 0.5 + 0.5*s;
pSafe = min(max(p, realmin), 1-eps);

nll = -sum(y.*log(pSafe) + (1-y).*log(1-pSafe));

dp = 0.5*s.*(1-s);
dEta = (p-y).*dp ./ (pSafe.*(1-pSafe));

gradient = zeros(nSessions+1,1);
gradient(1) = sum(dEta.*x);
gradient(2:end) = accumarray(sessionIndex, dEta, ...
  [nSessions 1], @sum, 0);
end
