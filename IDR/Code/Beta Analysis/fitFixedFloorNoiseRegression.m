function fit = fitFixedFloorNoiseRegression(xPref, xProbe, correct, varargin)
% fitFixedFloorNoiseRegression  Joint two-predictor logistic fit with 0.5 floor.
%
% Model:
%   P(correct) = 0.5 + 0.5*logistic(alpha + betaPref*xPref + betaProbe*xProbe)
%
% Predictors must be expressed in physical coherence-percentage-point units.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Display', 'off', @(x) ischar(x) || isstring(x));
addParameter(p, 'InitialBeta', [0.1 0.1], @(x) isnumeric(x) && numel(x) == 2 && all(isfinite(x)));
parse(p, varargin{:});

xPref = double(xPref(:));
xProbe = double(xProbe(:));
y = double(correct(:));

valid = isfinite(xPref) & isfinite(xProbe) & isfinite(y) & ismember(y, [0 1]);
xPref = xPref(valid);
xProbe = xProbe(valid);
y = y(valid);

nTrials = numel(y);
if nTrials < 3
  error('fitFixedFloorNoiseRegression:TooFewTrials', ...
    'At least three valid trials are required.');
end
nCorrect = sum(y == 1);
nError = sum(y == 0);
if nCorrect == 0 || nError == 0
  error('fitFixedFloorNoiseRegression:EmptyOutcomeClass', ...
    'The fit requires both correct and error trials.');
end

xCenter = [mean(xPref), mean(xProbe)];
X = [xPref - xCenter(1), xProbe - xCenter(2)];

pCorrect = mean(y);
q = min(max(2*pCorrect - 1, 1e-4), 1 - 1e-4);
alpha0 = log(q / (1-q));
theta0 = [double(p.Results.InitialBeta(:)); alpha0];

objective = @(theta) localNLLGradient(theta, X, y);
options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', char(p.Results.Display), ...
  'MaxFunctionEvaluations', 1e5, ...
  'MaxIterations', 2000, ...
  'OptimalityTolerance', 1e-9, ...
  'FunctionTolerance', 1e-12, ...
  'StepTolerance', 1e-12);

[thetaC, nll, exitflag, output, gradient, hessian] = ...
  fminunc(objective, theta0, options);

thetaC = full(double(thetaC));
gradient = full(double(gradient));
hessian = full(double(hessian));

beta = thetaC(1:2);
alphaCentered = thetaC(3);
alpha = alphaCentered - xCenter * beta;
theta = [beta; alpha];

covCentered = full(double(pinv(hessian)));
% [betaPref betaProbe alphaCentered] -> [betaPref betaProbe alpha]
J = eye(3);
J(3,1:2) = -xCenter;
covariance = full(double(J * covCentered * J'));
se = sqrt(max(0, diag(covariance)));

betaPref = beta(1);
betaProbe = beta(2);
ratio = betaProbe / betaPref;
ratioGradient = [-betaProbe / betaPref^2; 1 / betaPref];
betaCovariance = covariance(1:2, 1:2);
ratioVariance = ratioGradient' * betaCovariance * ratioGradient;
ratioSE = sqrt(max(0, ratioVariance));
ratioCI95 = ratio + 1.96 * ratioSE * [-1 1];

eta = alpha + betaPref*xPref + betaProbe*xProbe;
fittedProbability = localPsychometricProbability(eta);
logLikelihood = -nll;
nParameters = 3;

fit = struct();
fit.version = 1;
fit.model = ['P(correct)=0.5+0.5*logistic(' ...
  'alpha+betaPref*xPref+betaProbe*xProbe)'];
fit.chanceFloor = 0.5;
fit.nTrials = nTrials;
fit.nCorrect = nCorrect;
fit.nError = nError;
fit.fractionCorrect = mean(y);
fit.betaPref = betaPref;
fit.betaProbe = betaProbe;
fit.beta = beta;
fit.alpha = alpha;
fit.betaPrefSE = se(1);
fit.betaProbeSE = se(2);
fit.alphaSE = se(3);
fit.betaPrefCI95 = betaPref + 1.96*se(1)*[-1 1];
fit.betaProbeCI95 = betaProbe + 1.96*se(2)*[-1 1];
fit.alphaCI95 = alpha + 1.96*se(3)*[-1 1];
fit.betaRatio = ratio;
fit.betaRatioSE = ratioSE;
fit.betaRatioCI95 = ratioCI95;
fit.betaCovariance = betaCovariance;
fit.covariance = covariance;
fit.predictorCenter = xCenter;
fit.xPref = xPref;
fit.xProbe = xProbe;
fit.correct = logical(y);
fit.fittedProbability = fittedProbability;
fit.meanFittedProbability = mean(fittedProbability);
fit.logLikelihood = logLikelihood;
fit.AIC = 2*nParameters - 2*logLikelihood;
fit.BIC = log(nTrials)*nParameters - 2*logLikelihood;
fit.exitflag = exitflag;
fit.optimizerOutput = output;
fit.gradient = gradient;
fit.gradientInfNorm = norm(gradient, Inf);
fit.hessian = hessian;
fit.fitUsable = all(isfinite(theta)) && all(isfinite(se)) && ...
  isfinite(ratio) && isfinite(ratioSE);
end

function [nll, gradient] = localNLLGradient(theta, X, y)
beta = theta(1:2);
alpha = theta(3);
eta = alpha + X*beta;
[s, probability] = localLogisticAndProbability(eta);
pSafe = min(max(probability, realmin), 1-eps);
nll = -sum(y.*log(pSafe) + (1-y).*log(1-pSafe));
dp = 0.5*s.*(1-s);
dEta = (probability-y).*dp ./ (pSafe.*(1-pSafe));
gradient = [X' * dEta; sum(dEta)];
end

function [s, p] = localLogisticAndProbability(eta)
s = zeros(size(eta));
pos = eta >= 0;
s(pos) = 1 ./ (1 + exp(-eta(pos)));
e = exp(eta(~pos));
s(~pos) = e ./ (1 + e);
p = 0.5 + 0.5*s;
end

function p = localPsychometricProbability(eta)
[~, p] = localLogisticAndProbability(eta);
end
