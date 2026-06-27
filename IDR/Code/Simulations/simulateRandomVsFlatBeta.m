function sim = simulateRandomVsFlatBeta(nSessions, nTrialsPerSession, nSim, varargin)
% simulateRandomVsFlatBeta
% Compare across-session fitted betas for flat versus random 3-direction readout.
%
% The two models are matched for overall performance.
%
% FLAT:
%   All three directions contribute equally on every trial.
%   One common effective operating point applies on every trial.
%
% RANDOM:
%   One direction is monitored per trial, selected uniformly.
%   When the monitored direction matches the stepped direction, the signal
%   is strong; otherwise it is attenuated by r120.
%
% Both models use the pooled psychometric slope and fixed 0.5 floor:
%
%   P(correct) = 0.5 + 0.5*logistic(alpha + beta*x)
%
% The session-level beta is estimated from simulated preferred-noise trials
% using the same single-session fitting model used elsewhere.
%
% Usage:
%   sim = simulateRandomVsFlatBeta(20, 1000, 1000);
%
% Optional name/value arguments:
%   'TargetPerformance'   default 0.75
%   'R120'                default 0.20
%   'Seed'                default 1
%   'BetaTrue'            default pooled psychometric beta
%   'SessionIntercepts'   default pooled psychometric session intercepts
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/
%       RandomVsFlatBetaSimulation.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>

P = inputParser;
addParameter(P, 'TargetPerformance', 0.75, ...
  @(x) isnumeric(x) && isscalar(x) && x > 0.5 && x < 1);
addParameter(P, 'R120', 0.20, ...
  @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(P, 'Seed', 1, ...
  @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(P, 'BetaTrue', [], ...
  @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x)));
addParameter(P, 'SessionIntercepts', [], ...
  @(x) isempty(x) || isnumeric(x));
parse(P, varargin{:});
R = P.Results;

validateattributes(nSessions, {'numeric'}, ...
  {'scalar','integer','positive'});
validateattributes(nTrialsPerSession, {'numeric'}, ...
  {'scalar','integer','positive'});
validateattributes(nSim, {'numeric'}, ...
  {'scalar','integer','positive'});

rng(R.Seed);

baseFolder = domainFolder(mfilename('fullpath'));
acrossFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions');

psychPath = fullfile(acrossFolder, 'BetaPsychometricFit.mat');
S = load(psychPath, 'psychometricFit');
F = S.psychometricFit;

if isempty(R.BetaTrue)
  betaTrue = double(F.beta);
else
  betaTrue = double(R.BetaTrue);
end

if isempty(R.SessionIntercepts)
  interceptPool = double(F.alpha(:));
else
  interceptPool = double(R.SessionIntercepts(:));
end

targetPerf = R.TargetPerformance;
r120 = R.R120;

% Use the mean session intercept to define the physical step coherence
% needed for each model. Individual simulated sessions draw intercepts from
% the empirical pool.
alphaRef = mean(interceptPool);

% Flat model: one common effective signal.
flatEffectiveSignal = invertPsychometric(targetPerf, alphaRef, betaTrue);

% Under equal weighting of three directions, physical step coherence is
% three times the effective signal.
flatPhysicalStep = 3 * flatEffectiveSignal;

% Random model:
%   state 1 probability = 1/3, effective signal = physical step
%   states 2-3 total probability = 2/3, effective signal = r120*physical step
randomPhysicalStep = solveRandomStep(targetPerf, alphaRef, betaTrue, r120);

randomSignals = [randomPhysicalStep, r120*randomPhysicalStep];
randomStateProb = [1/3, 2/3];
randomStatePerf = psychometricProbability( ...
  alphaRef + betaTrue*randomSignals);
randomAchievedPerf = sum(randomStateProb .* randomStatePerf);

flatStatePerf = psychometricProbability( ...
  alphaRef + betaTrue*flatEffectiveSignal);

% Simulate session fits.
flatSessionBetas = nan(nSim, nSessions);
randomSessionBetas = nan(nSim, nSessions);
flatMeanBeta = nan(nSim,1);
randomMeanBeta = nan(nSim,1);
meanDifference = nan(nSim,1);

for iSim = 1:nSim
  if iSim == 1 || mod(iSim,100) == 0
    fprintf('Simulation %d of %d\n', iSim, nSim);
  end

  alphaDraw = interceptPool(randi(numel(interceptPool), nSessions, 1));

  for s = 1:nSessions
    % ---- Flat ----
    signalFlat = repmat(flatEffectiveSignal, nTrialsPerSession, 1);
    noiseFlat = randn(nTrialsPerSession,1);
    noiseFlat = noiseFlat / std(noiseFlat);

    etaFlat = alphaDraw(s) + betaTrue*(signalFlat + noiseFlat);
    yFlat = double(rand(nTrialsPerSession,1) < ...
      psychometricProbability(etaFlat));

    [flatSessionBetas(iSim,s), ~] = fitOneSession(noiseFlat, yFlat);

    % ---- Random ----
    selectedHigh = rand(nTrialsPerSession,1) < 1/3;
    signalRandom = r120*randomPhysicalStep * ...
      ones(nTrialsPerSession,1);
    signalRandom(selectedHigh) = randomPhysicalStep;

    noiseRandom = randn(nTrialsPerSession,1);
    noiseRandom = noiseRandom / std(noiseRandom);

    % Noise is behaviorally effective only when its direction is the
    % monitored direction. This creates the random-state mixture.
    noiseGain = r120 * ones(nTrialsPerSession,1);
    noiseGain(selectedHigh) = 1;

    etaRandom = alphaDraw(s) + ...
      betaTrue*(signalRandom + noiseGain.*noiseRandom);

    yRandom = double(rand(nTrialsPerSession,1) < ...
      psychometricProbability(etaRandom));

    [randomSessionBetas(iSim,s), ~] = ...
      fitOneSession(noiseRandom, yRandom);
  end

  flatMeanBeta(iSim) = mean(flatSessionBetas(iSim,:), 'omitnan');
  randomMeanBeta(iSim) = mean(randomSessionBetas(iSim,:), 'omitnan');
  meanDifference(iSim) = ...
    randomMeanBeta(iSim) - flatMeanBeta(iSim);
end

sim = struct();
sim.version = 1;
sim.nSessions = nSessions;
sim.nTrialsPerSession = nTrialsPerSession;
sim.nSim = nSim;
sim.seed = R.Seed;
sim.betaTrue = betaTrue;
sim.targetPerformance = targetPerf;
sim.r120 = r120;

sim.flat = struct();
sim.flat.physicalStepCohPC = flatPhysicalStep;
sim.flat.effectiveSignalCohPC = flatEffectiveSignal;
sim.flat.stateProbabilities = 1;
sim.flat.statePerformances = flatStatePerf;
sim.flat.achievedPerformance = flatStatePerf;
sim.flat.sessionBetas = flatSessionBetas;
sim.flat.meanBetas = flatMeanBeta;
sim.flat.meanBeta = mean(flatMeanBeta, 'omitnan');
sim.flat.meanBetaCI95 = prctile(flatMeanBeta, [2.5 97.5]);

sim.random = struct();
sim.random.physicalStepCohPC = randomPhysicalStep;
sim.random.effectiveSignalCohPC = randomSignals;
sim.random.stateProbabilities = randomStateProb;
sim.random.statePerformances = randomStatePerf;
sim.random.achievedPerformance = randomAchievedPerf;
sim.random.sessionBetas = randomSessionBetas;
sim.random.meanBetas = randomMeanBeta;
sim.random.meanBeta = mean(randomMeanBeta, 'omitnan');
sim.random.meanBetaCI95 = prctile(randomMeanBeta, [2.5 97.5]);

sim.comparison = struct();
sim.comparison.meanDifference = mean(meanDifference, 'omitnan');
sim.comparison.differenceCI95 = prctile(meanDifference, [2.5 97.5]);
sim.comparison.probRandomBelowFlat = mean(meanDifference < 0, 'omitnan');
sim.comparison.meanDifferenceSamples = meanDifference;

outputPath = fullfile(acrossFolder, ...
  'RandomVsFlatBetaSimulation.mat');
save(outputPath, 'sim', '-v7.3');

fprintf('\nTarget overall performance: %.4f\n', targetPerf);

fprintf('\nFLAT READOUT\n');
fprintf('  Physical step coherence: %.4f%%\n', flatPhysicalStep);
fprintf('  Effective signal coherence: %.4f%%\n', flatEffectiveSignal);
fprintf('  Operating-point performance: %.4f\n', flatStatePerf);
fprintf('  Mean fitted beta: %.6g (95%% CI %.6g to %.6g)\n', ...
  sim.flat.meanBeta, sim.flat.meanBetaCI95(1), ...
  sim.flat.meanBetaCI95(2));

fprintf('\nRANDOM READOUT\n');
fprintf('  Physical step coherence: %.4f%%\n', randomPhysicalStep);
fprintf('  High-state signal: %.4f%%, probability %.3f, performance %.4f\n', ...
  randomSignals(1), randomStateProb(1), randomStatePerf(1));
fprintf('  Low-state signal: %.4f%%, probability %.3f, performance %.4f\n', ...
  randomSignals(2), randomStateProb(2), randomStatePerf(2));
fprintf('  Achieved overall performance: %.4f\n', randomAchievedPerf);
fprintf('  Mean fitted beta: %.6g (95%% CI %.6g to %.6g)\n', ...
  sim.random.meanBeta, sim.random.meanBetaCI95(1), ...
  sim.random.meanBetaCI95(2));

fprintf('\nCOMPARISON\n');
fprintf('  Random - flat mean beta: %.6g\n', ...
  sim.comparison.meanDifference);
fprintf('  95%% interval: %.6g to %.6g\n', ...
  sim.comparison.differenceCI95(1), ...
  sim.comparison.differenceCI95(2));
fprintf('  P(random mean beta < flat mean beta): %.4f\n', ...
  sim.comparison.probRandomBelowFlat);
fprintf('Saved %s\n', outputPath);

plotDistributions(sim);
end

% -------------------------------------------------------------------------
function step = solveRandomStep(targetPerf, alpha, beta, r120)

objective = @(step) ...
  (1/3)*psychometricProbability(alpha + beta*step) + ...
  (2/3)*psychometricProbability(alpha + beta*r120*step) - ...
  targetPerf;

lo = 0;
hi = 100;

while objective(hi) < 0
  hi = hi*2;
end

step = fzero(objective, [lo hi]);
end

% -------------------------------------------------------------------------
function x = invertPsychometric(targetPerf, alpha, beta)

q = 2*targetPerf - 1;
eta = log(q/(1-q));
x = (eta-alpha)/beta;
end

% -------------------------------------------------------------------------
function [betaHat, betaSE] = fitOneSession(x, y)

x = double(x(:));
y = double(y(:));

if all(y == 1) || all(y == 0)
  betaHat = NaN;
  betaSE = NaN;
  return;
end

xCenter = mean(x);
xC = x-xCenter;

pObs = mean(y);
q = min(max(2*pObs-1,1e-4),1-1e-4);
theta0 = [0.1; log(q/(1-q))];

objective = @(theta) nllGrad(theta,xC,y);

options = optimoptions('fminunc', ...
  'Algorithm','trust-region', ...
  'SpecifyObjectiveGradient',true, ...
  'Display','off', ...
  'MaxFunctionEvaluations',2e4, ...
  'MaxIterations',500, ...
  'FunctionTolerance',1e-10, ...
  'OptimalityTolerance',1e-7, ...
  'StepTolerance',1e-10);

[theta,~,exitflag,~,~,hessian] = ...
  fminunc(objective,theta0,options);

if exitflag <= 0
  betaHat = NaN;
  betaSE = NaN;
  return;
end

hessian = full(double(hessian));
covTheta = full(inv(hessian));

betaHat = full(double(theta(1)));
betaSE = sqrt(covTheta(1,1));
end

% -------------------------------------------------------------------------
function [nll,gradient] = nllGrad(theta,x,y)

beta = theta(1);
alpha = theta(2);

eta = alpha+beta*x;
[s,p] = logisticAndPsychometricProbability(eta);

pSafe = min(max(p,realmin),1-eps);
nll = -sum(y.*log(pSafe)+(1-y).*log(1-pSafe));

dp = 0.5*s.*(1-s);
dEta = (p-y).*dp./(pSafe.*(1-pSafe));

gradient = [sum(dEta.*x); sum(dEta)];
end

% -------------------------------------------------------------------------
function p = psychometricProbability(eta)
[~,p] = logisticAndPsychometricProbability(eta);
end

% -------------------------------------------------------------------------
function [s,p] = logisticAndPsychometricProbability(eta)

s = zeros(size(eta));
positive = eta >= 0;
s(positive) = 1./(1+exp(-eta(positive)));

z = exp(eta(~positive));
s(~positive) = z./(1+z);

p = 0.5+0.5*s;
end

% -------------------------------------------------------------------------
function plotDistributions(sim)

figure;
histogram(sim.flat.meanBetas);
hold on
histogram(sim.random.meanBetas);
xlabel('Across-session mean fitted beta');
ylabel('Number of simulations');
legend('Flat','Random');
title(sprintf('%d sessions, %d trials/session', ...
  sim.nSessions, sim.nTrialsPerSession));
box off
end
