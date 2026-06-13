function precisionSim = simulateBetaPrecisionByTrialCount(trialCounts, nSim, seed)
% simulateBetaPrecisionByTrialCount
% Estimate single-session beta precision as a function of trial count.
%
% For each target trial count and simulation:
%   - choose one of the observed sessions
%   - sample that session's actual effectiveNoisePC values with replacement
%   - simulate outcomes using the pooled common beta and that session's
%     pooled-fit intercept
%   - fit one session beta
%   - compare the fitted beta with the pooled common beta
%
% This isolates finite-sample precision while preserving realistic noise
% distributions and session operating points.
%
% Usage:
%   precisionSim = simulateBetaPrecisionByTrialCount();
%   precisionSim = simulateBetaPrecisionByTrialCount([250 500 750 1000],1000,1);
%
% Reads:
%   BetaNoiseRegressionFit.mat
%   SessionData/*.mat
%
% Saves:
%   BetaPrecisionByTrialCount.mat

cleanupObj = initProjectPath(); %#ok<NASGU>

if nargin < 1 || isempty(trialCounts)
  trialCounts = [250 500 750 1000];
end
if nargin < 2 || isempty(nSim)
  nSim = 1000;
end
if nargin < 3 || isempty(seed)
  seed = 1;
end

trialCounts = double(trialCounts(:)');
rng(seed);

baseFolder = folderPath();
sessionFolder = fullfile(baseFolder,'Data','FullSessions', ...
  'BetaAnalysis','SessionData');
acrossFolder = fullfile(baseFolder,'Data','FullSessions', ...
  'BetaAnalysis','AcrossSessions');

P = load(fullfile(acrossFolder,'BetaNoiseRegressionFit.mat'),'noiseFit');
F = P.noiseFit;

betaTrue = double(F.betaNoise);
fileNames = string(F.sessionFileNames(:));
alpha = double(F.alpha(:));
nSessions = numel(fileNames);

sessionX = cell(nSessions,1);

for s = 1:nSessions
  D = load(fullfile(sessionFolder,fileNames(s)),'sessionNoise');
  N = D.sessionNoise;
  use = logical(N.hasPreferredNoise(:));
  x = double(N.trialAnalysis.effectiveNoisePC(use));
  x = x(:);

  if isempty(x) || any(~isfinite(x))
    error('simulateBetaPrecisionByTrialCount:BadSessionData', ...
      '%s has invalid preferred-noise predictor values.', fileNames(s));
  end

  sessionX{s} = x;
end

nCounts = numel(trialCounts);

betaHat = nan(nSim,nCounts);
betaSE = nan(nSim,nCounts);
covered = false(nSim,nCounts);
sourceSession = nan(nSim,nCounts);
nError = nan(nSim,nCounts);
fractionCorrect = nan(nSim,nCounts);

for c = 1:nCounts
  nTrials = trialCounts(c);

  fprintf('\nTrial count %d\n', nTrials);

  for b = 1:nSim
    if b == 1 || mod(b,100) == 0
      fprintf('  simulation %d of %d\n', b, nSim);
    end

    s = randi(nSessions);
    xSource = sessionX{s};

    idx = randi(numel(xSource),nTrials,1);
    x = xSource(idx);

    eta = alpha(s) + betaTrue*x;
    p = psychometricProbability(eta);
    y = double(rand(nTrials,1) < p);

    % All-correct or all-error samples cannot support a slope fit.
    if all(y == 1) || all(y == 0)
      continue;
    end

    [bHat,seHat,exitflag] = fitOneSession(x,y);

    if exitflag <= 0 || ~isfinite(bHat) || ~isfinite(seHat) || seHat <= 0
      continue;
    end

    betaHat(b,c) = bHat;
    betaSE(b,c) = seHat;
    covered(b,c) = betaTrue >= bHat-1.96*seHat && ...
                   betaTrue <= bHat+1.96*seHat;
    sourceSession(b,c) = s;
    nError(b,c) = sum(y == 0);
    fractionCorrect(b,c) = mean(y);
  end
end

summaryTable = table('Size',[nCounts 12], ...
  'VariableTypes',repmat({'double'},1,12), ...
  'VariableNames',{'nTrials','nValid','meanBeta','medianBeta','bias', ...
  'sdBeta','rmse','medianAbsError','p80AbsError','p90AbsError', ...
  'coverage95','medianSE'});

for c = 1:nCounts
  valid = isfinite(betaHat(:,c));
  err = betaHat(valid,c)-betaTrue;
  absErr = abs(err);

  summaryTable.nTrials(c) = trialCounts(c);
  summaryTable.nValid(c) = sum(valid);
  summaryTable.meanBeta(c) = mean(betaHat(valid,c));
  summaryTable.medianBeta(c) = median(betaHat(valid,c));
  summaryTable.bias(c) = mean(err);
  summaryTable.sdBeta(c) = std(betaHat(valid,c));
  summaryTable.rmse(c) = sqrt(mean(err.^2));
  summaryTable.medianAbsError(c) = median(absErr);
  summaryTable.p80AbsError(c) = prctile(absErr,80);
  summaryTable.p90AbsError(c) = prctile(absErr,90);
  summaryTable.coverage95(c) = mean(covered(valid,c));
  summaryTable.medianSE(c) = median(betaSE(valid,c));
end

precisionSim = struct();
precisionSim.version = 1;
precisionSim.betaTrue = betaTrue;
precisionSim.trialCounts = trialCounts;
precisionSim.nSim = nSim;
precisionSim.seed = seed;
precisionSim.betaHat = betaHat;
precisionSim.betaSE = betaSE;
precisionSim.covered = covered;
precisionSim.sourceSession = sourceSession;
precisionSim.nError = nError;
precisionSim.fractionCorrect = fractionCorrect;
precisionSim.summaryTable = summaryTable;
precisionSim.createdBy = mfilename;
precisionSim.createdDate = datetime('now');

outputPath = fullfile(acrossFolder,'BetaPrecisionByTrialCount.mat');
save(outputPath,'precisionSim','-v7.3');

disp(summaryTable);
fprintf('True pooled noise beta: %.6g\n', betaTrue);
fprintf('Saved %s\n', outputPath);

plotPrecision(summaryTable,betaTrue);
end

% -------------------------------------------------------------------------
function [beta,betaSE,exitflag] = fitOneSession(x,y)

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

[theta,~,exitflag,~,~,hessian] = fminunc(objective,theta0,options);

if exitflag <= 0
  beta = NaN;
  betaSE = NaN;
  return;
end

hessian = full(double(hessian));
covTheta = full(inv(hessian));

beta = full(double(theta(1)));
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
function plotPrecision(T,betaTrue)

figure;
plot(T.nTrials,T.rmse,'ko-','LineWidth',1.5,'MarkerFaceColor','w');
xlabel('Preferred-noise trials per session');
ylabel('RMSE of fitted beta');
title(sprintf('Single-session beta precision; true beta %.4g',betaTrue));
box off

figure;
plot(T.nTrials,T.medianAbsError,'ko-','LineWidth',1.5,'MarkerFaceColor','w');
hold on
plot(T.nTrials,T.p80AbsError,'k--');
plot(T.nTrials,T.p90AbsError,'k:');
xlabel('Preferred-noise trials per session');
ylabel('Absolute beta error');
legend('Median','80th percentile','90th percentile','Location','best');
title('Expected single-session beta error');
box off
end
