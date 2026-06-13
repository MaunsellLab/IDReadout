function sim = simulateRandomVsFlatBetaEmpirical(nSessions, nTrialsPerSession, nSim, varargin)
% simulateRandomVsFlatBetaEmpirical
% Definitive random-versus-flat simulation using empirical effective-noise
% distributions and one fitted beta for the mean of three noise streams.
%
% Each simulated trial contains three independent noise streams sampled
% from one real session's empirical effectiveNoisePC distribution.
%
% Analysis predictor:
%   meanNoisePC = mean([noise1 noise2 noise3], 2)
%
% FLAT READOUT
%   Decision evidence:
%       physicalStep/3 + meanNoisePC
%
% RANDOM READOUT
%   One latent monitored direction is selected uniformly each trial.
%   The stepped direction is also selected uniformly and independently.
%   The readout weights are peaked:
%       weight = 1 for the monitored direction
%       weight = r120 for either other direction
%
%   Decision evidence:
%       signalWeight*physicalStep + ...
%       sum(readoutWeights .* noiseStreams)
%
% Both models use:
%   P(correct) = 0.5 + 0.5*logistic(alpha_session + beta*evidence)
%
% The physical step is calibrated separately for every empirical session
% and model so that expected overall performance, including coherence noise,
% equals TargetPerformance.
%
% Each simulated session is analyzed exactly as planned:
%   Correct ~ alpha + betaMean*meanNoisePC
%
% Usage:
%   sim = simulateRandomVsFlatBetaEmpirical(20, 1000, 1000);
%
% Optional parameters:
%   'TargetPerformance'      default 0.75
%   'R120'                   default 0.20
%   'Seed'                   default 1
%   'CalibrationSamples'     default 50000
%
% Saves:
%   Data/FullSessions/BetaAnalysis/AcrossSessions/
%       RandomVsFlatBetaEmpiricalSimulation.mat

cleanupObj = initProjectPath(); %#ok<NASGU>

P = inputParser;
addParameter(P, 'TargetPerformance', 0.75, ...
  @(x) isnumeric(x) && isscalar(x) && x > 0.5 && x < 1);
addParameter(P, 'R120', 0.20, ...
  @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
addParameter(P, 'Seed', 1, ...
  @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(P, 'CalibrationSamples', 50000, ...
  @(x) isnumeric(x) && isscalar(x) && x >= 1000);
parse(P, varargin{:});
R = P.Results;

validateattributes(nSessions, {'numeric'}, ...
  {'scalar','integer','positive'});
validateattributes(nTrialsPerSession, {'numeric'}, ...
  {'scalar','integer','positive'});
validateattributes(nSim, {'numeric'}, ...
  {'scalar','integer','positive'});

rng(R.Seed);

baseFolder = folderPath();
sessionFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'SessionData');
acrossFolder = fullfile(baseFolder, 'Data', 'FullSessions', ...
  'BetaAnalysis', 'AcrossSessions');

Psy = load(fullfile(acrossFolder, 'BetaPsychometricFit.mat'), ...
  'psychometricFit');
F = Psy.psychometricFit;

betaTrue = double(F.beta);
fileNames = string(F.sessionFileNames(:));
alpha = double(F.alpha(:));
nSourceSessions = numel(fileNames);

noisePool = cell(nSourceSessions,1);

for s = 1:nSourceSessions
  D = load(fullfile(sessionFolder, fileNames(s)), 'sessionNoise');
  N = D.sessionNoise;
  use = logical(N.hasPreferredNoise(:));
  x = double(N.trialAnalysis.effectiveNoisePC(use));
  x = x(:);

  if isempty(x) || any(~isfinite(x))
    error('simulateRandomVsFlatBetaEmpirical:BadNoisePool', ...
      '%s has invalid effective-noise values.', fileNames(s));
  end

  noisePool{s} = x;
end

target = R.TargetPerformance;
r120 = R.R120;
nCal = R.CalibrationSamples;

% Precompute one calibrated step per empirical source session and model.
flatStep = nan(nSourceSessions,1);
randomStep = nan(nSourceSessions,1);

flatOperatingSignal = nan(nSourceSessions,1);
flatOperatingPerformanceNoNoise = nan(nSourceSessions,1);
flatAchievedPerformance = nan(nSourceSessions,1);

randomHighSignal = nan(nSourceSessions,1);
randomLowSignal = nan(nSourceSessions,1);
randomHighPerformanceNoNoise = nan(nSourceSessions,1);
randomLowPerformanceNoNoise = nan(nSourceSessions,1);
randomAchievedPerformance = nan(nSourceSessions,1);

for s = 1:nSourceSessions
  xPool = noisePool{s};

  noise3 = sampleThreeStreams(xPool, nCal);
  meanNoise = mean(noise3,2);

  % Fixed latent states for calibration make the objective deterministic.
  stepIndex = randi(3,nCal,1);
  monitorIndex = randi(3,nCal,1);
  randomNoiseEvidence = randomEvidenceNoise( ...
    noise3, monitorIndex, r120);
  signalGain = r120*ones(nCal,1);
  signalGain(monitorIndex == stepIndex) = 1;

  flatObjective = @(step) mean(psychometricProbability( ...
    alpha(s) + betaTrue*(step/3 + meanNoise))) - target;

  randomObjective = @(step) mean(psychometricProbability( ...
    alpha(s) + betaTrue*(signalGain*step + randomNoiseEvidence))) ...
    - target;

  flatStep(s) = solvePositiveRoot(flatObjective);
  randomStep(s) = solvePositiveRoot(randomObjective);

  flatOperatingSignal(s) = flatStep(s)/3;
  flatOperatingPerformanceNoNoise(s) = psychometricProbability( ...
    alpha(s) + betaTrue*flatOperatingSignal(s));
  flatAchievedPerformance(s) = flatObjective(flatStep(s)) + target;

  randomHighSignal(s) = randomStep(s);
  randomLowSignal(s) = r120*randomStep(s);
  randomHighPerformanceNoNoise(s) = psychometricProbability( ...
    alpha(s) + betaTrue*randomHighSignal(s));
  randomLowPerformanceNoNoise(s) = psychometricProbability( ...
    alpha(s) + betaTrue*randomLowSignal(s));
  randomAchievedPerformance(s) = ...
    randomObjective(randomStep(s)) + target;
end

flatSessionBeta = nan(nSim,nSessions);
randomSessionBeta = nan(nSim,nSessions);
flatMeanBeta = nan(nSim,1);
randomMeanBeta = nan(nSim,1);
difference = nan(nSim,1);
sourceSession = nan(nSim,nSessions);

for iSim = 1:nSim
  if iSim == 1 || mod(iSim,50) == 0
    fprintf('Simulation %d of %d\n', iSim, nSim);
  end

  source = randi(nSourceSessions,nSessions,1);
  sourceSession(iSim,:) = source;

  for j = 1:nSessions
    s = source(j);
    xPool = noisePool{s};

    % ----- Flat -----
    noise3 = sampleThreeStreams(xPool,nTrialsPerSession);
    meanNoise = mean(noise3,2);

    evidenceFlat = flatStep(s)/3 + meanNoise;
    pFlat = psychometricProbability( ...
      alpha(s) + betaTrue*evidenceFlat);
    yFlat = double(rand(nTrialsPerSession,1) < pFlat);

    flatSessionBeta(iSim,j) = fitOneSession(meanNoise,yFlat);

    % ----- Random -----
    noise3 = sampleThreeStreams(xPool,nTrialsPerSession);
    meanNoise = mean(noise3,2);

    stepIndex = randi(3,nTrialsPerSession,1);
    monitorIndex = randi(3,nTrialsPerSession,1);

    signalGain = r120*ones(nTrialsPerSession,1);
    signalGain(monitorIndex == stepIndex) = 1;

    noiseEvidence = randomEvidenceNoise( ...
      noise3, monitorIndex, r120);

    evidenceRandom = signalGain*randomStep(s) + noiseEvidence;
    pRandom = psychometricProbability( ...
      alpha(s) + betaTrue*evidenceRandom);
    yRandom = double(rand(nTrialsPerSession,1) < pRandom);

    randomSessionBeta(iSim,j) = ...
      fitOneSession(meanNoise,yRandom);
  end

  flatMeanBeta(iSim) = mean(flatSessionBeta(iSim,:),'omitnan');
  randomMeanBeta(iSim) = mean(randomSessionBeta(iSim,:),'omitnan');
  difference(iSim) = randomMeanBeta(iSim)-flatMeanBeta(iSim);
end

sim = struct();
sim.version = 1;
sim.modelDescription = ...
  ['three independent empirical effective-noise streams; ' ...
   'analysis regressor is their mean'];
sim.nSessions = nSessions;
sim.nTrialsPerSession = nTrialsPerSession;
sim.nSim = nSim;
sim.seed = R.Seed;
sim.targetPerformance = target;
sim.r120 = r120;
sim.betaTrue = betaTrue;
sim.sourceSessionFileNames = cellstr(fileNames);

sim.operatingPoints = struct();

sim.operatingPoints.flat = table( ...
  cellstr(fileNames), flatStep, flatOperatingSignal, ...
  flatOperatingPerformanceNoNoise, flatAchievedPerformance, ...
  'VariableNames', {'fileName','physicalStepCohPC', ...
  'effectiveSignalCohPC','noNoisePerformance','achievedPerformance'});

sim.operatingPoints.random = table( ...
  cellstr(fileNames), randomStep, randomHighSignal, randomLowSignal, ...
  randomHighPerformanceNoNoise, randomLowPerformanceNoNoise, ...
  randomAchievedPerformance, ...
  'VariableNames', {'fileName','physicalStepCohPC', ...
  'highSignalCohPC','lowSignalCohPC', ...
  'highNoNoisePerformance','lowNoNoisePerformance', ...
  'achievedPerformance'});

sim.flat.sessionBetas = flatSessionBeta;
sim.flat.meanBetas = flatMeanBeta;
sim.flat.meanBeta = mean(flatMeanBeta,'omitnan');
sim.flat.meanBetaCI95 = prctile(flatMeanBeta,[2.5 97.5]);

sim.random.sessionBetas = randomSessionBeta;
sim.random.meanBetas = randomMeanBeta;
sim.random.meanBeta = mean(randomMeanBeta,'omitnan');
sim.random.meanBetaCI95 = prctile(randomMeanBeta,[2.5 97.5]);

sim.comparison.differenceSamples = difference;
sim.comparison.meanDifference = mean(difference,'omitnan');
sim.comparison.differenceCI95 = prctile(difference,[2.5 97.5]);
sim.comparison.probRandomBelowFlat = mean(difference < 0,'omitnan');
sim.sourceSession = sourceSession;
sim.createdBy = mfilename;
sim.createdDate = datetime('now');

outputPath = fullfile(acrossFolder, ...
  'RandomVsFlatBetaEmpiricalSimulation.mat');
save(outputPath,'sim','-v7.3');

fprintf('\nTarget overall performance: %.4f\n',target);
fprintf('True psychometric beta: %.6g\n',betaTrue);
fprintf('r120: %.3f\n',r120);

fprintf('\nFLAT OPERATING POINTS ACROSS SOURCE SESSIONS\n');
printRange('Physical step coherence',flatStep,'%%');
printRange('Effective signal coherence',flatOperatingSignal,'%%');
printRange('No-noise performance',flatOperatingPerformanceNoNoise,'');
printRange('Achieved performance',flatAchievedPerformance,'');

fprintf('\nRANDOM OPERATING POINTS ACROSS SOURCE SESSIONS\n');
printRange('Physical step coherence',randomStep,'%%');
printRange('High-state signal',randomHighSignal,'%%');
printRange('Low-state signal',randomLowSignal,'%%');
printRange('High-state no-noise performance', ...
  randomHighPerformanceNoNoise,'');
printRange('Low-state no-noise performance', ...
  randomLowPerformanceNoNoise,'');
printRange('Achieved performance',randomAchievedPerformance,'');

fprintf('\nFITTED ACROSS-SESSION BETAS\n');
fprintf('Flat mean beta: %.6g (95%% interval %.6g to %.6g)\n', ...
  sim.flat.meanBeta,sim.flat.meanBetaCI95(1),sim.flat.meanBetaCI95(2));
fprintf('Random mean beta: %.6g (95%% interval %.6g to %.6g)\n', ...
  sim.random.meanBeta,sim.random.meanBetaCI95(1),sim.random.meanBetaCI95(2));
fprintf('Random - flat: %.6g (95%% interval %.6g to %.6g)\n', ...
  sim.comparison.meanDifference, ...
  sim.comparison.differenceCI95(1), ...
  sim.comparison.differenceCI95(2));
fprintf('P(random < flat): %.4f\n', ...
  sim.comparison.probRandomBelowFlat);
fprintf('Saved %s\n',outputPath);

plotSimulation(sim);
end

% -------------------------------------------------------------------------
function X = sampleThreeStreams(pool,n)
idx = randi(numel(pool),n,3);
X = reshape(pool(idx),n,3);
end

% -------------------------------------------------------------------------
function e = randomEvidenceNoise(noise3,monitorIndex,r120)
n = size(noise3,1);
weights = r120*ones(n,3);
linear = sub2ind([n 3],(1:n)',monitorIndex);
weights(linear) = 1;
e = sum(weights.*noise3,2);
end

% -------------------------------------------------------------------------
function root = solvePositiveRoot(fun)
lo = 0;
hi = 50;
while fun(hi) < 0
  hi = hi*2;
  if hi > 10000
    error('simulateRandomVsFlatBetaEmpirical:NoRoot', ...
      'Could not bracket target-performance step.');
  end
end
root = fzero(fun,[lo hi]);
end

% -------------------------------------------------------------------------
function betaHat = fitOneSession(x,y)
x = double(x(:));
y = double(y(:));

if all(y==1) || all(y==0)
  betaHat = NaN;
  return;
end

xCenter = mean(x);
xC = x-xCenter;
pObs = mean(y);
q = min(max(2*pObs-1,1e-4),1-1e-4);
theta0 = [0.1;log(q/(1-q))];

objective = @(theta)nllGrad(theta,xC,y);
options = optimoptions('fminunc', ...
  'Algorithm','trust-region', ...
  'SpecifyObjectiveGradient',true, ...
  'Display','off', ...
  'MaxFunctionEvaluations',2e4, ...
  'MaxIterations',500, ...
  'FunctionTolerance',1e-10, ...
  'OptimalityTolerance',1e-7, ...
  'StepTolerance',1e-10);

[theta,~,exitflag] = fminunc(objective,theta0,options);

if exitflag <= 0
  betaHat = NaN;
else
  betaHat = full(double(theta(1)));
end
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
gradient = [sum(dEta.*x);sum(dEta)];
end

% -------------------------------------------------------------------------
function p = psychometricProbability(eta)
[~,p] = logisticAndPsychometricProbability(eta);
end

% -------------------------------------------------------------------------
function [s,p] = logisticAndPsychometricProbability(eta)
s = zeros(size(eta));
positive = eta>=0;
s(positive) = 1./(1+exp(-eta(positive)));
z = exp(eta(~positive));
s(~positive) = z./(1+z);
p = 0.5+0.5*s;
end

% -------------------------------------------------------------------------
function printRange(label,x,suffix)
fprintf('  %s: median %.4f%s, range %.4f to %.4f%s\n', ...
  label,median(x),suffix,min(x),max(x),suffix);
end

% -------------------------------------------------------------------------
function plotSimulation(sim)
figure;
histogram(sim.flat.meanBetas);
hold on
histogram(sim.random.meanBetas);
xlabel('Across-session mean fitted \beta_{mean}');
ylabel('Number of simulations');
legend('Flat','Random','Location','best');
title(sprintf('%d sessions, %d trials/session', ...
  sim.nSessions,sim.nTrialsPerSession));
box off
end
