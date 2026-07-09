function fit = fitSharedProbeScale(sessionData, varargin)
% fitSharedProbeScale  Estimate one probe/pref scale shared across sessions.
%
% For trial i in session s:
%   P(correct) = 0.5 + 0.5*logistic(alpha_s + ...
%                  betaPref_s*(xPref + scale*xProbe))
%
% sessionData is a cell array. Each cell must contain vectors:
%   xPref, xProbe, correct
%
% The session-specific preferred sensitivities are nuisance parameters.
% The shared scale is the normalized probe sensitivity of interest.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Display', 'off', @(x) ischar(x) || isstring(x));
addParameter(p, 'InitialScale', 0.3, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'InitialBetaPref', 0.1, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
parse(p, varargin{:});

sessionData = sessionData(:);
nSessions = numel(sessionData);
if nSessions < 1
  error('fitSharedProbeScale:NoSessions', 'At least one session is required.');
end

xPref = cell(nSessions,1);
xProbe = cell(nSessions,1);
y = cell(nSessions,1);
nTrialsBySession = zeros(nSessions,1);
nCorrectBySession = zeros(nSessions,1);

for s = 1:nSessions
  D = sessionData{s};
  required = {'xPref','xProbe','correct'};
  missing = required(~isfield(D, required));
  if ~isempty(missing)
    error('fitSharedProbeScale:MissingFields', ...
      'Session %d is missing: %s', s, strjoin(missing, ', '));
  end

  xp = double(D.xPref(:));
  xq = double(D.xProbe(:));
  yy = double(D.correct(:));
  valid = isfinite(xp) & isfinite(xq) & isfinite(yy) & ismember(yy,[0 1]);
  xp = xp(valid); xq = xq(valid); yy = yy(valid);

  if numel(yy) < 3 || all(yy == 0) || all(yy == 1)
    error('fitSharedProbeScale:UnusableSession', 'Session %d lacks enough trials or both outcome classes.', s);
  end

  xPref{s} = xp;
  xProbe{s} = xq;
  y{s} = yy;
  nTrialsBySession(s) = numel(yy);
  nCorrectBySession(s) = sum(yy);
end

% theta = [scale; betaPref_1..S; alpha_1..S]
scale0 = double(p.Results.InitialScale);
beta0 = repmat(double(p.Results.InitialBetaPref), nSessions, 1);
alpha0 = nan(nSessions,1);
for s = 1:nSessions
  pc = mean(y{s});
  q = min(max(2*pc - 1, 1e-4), 1 - 1e-4);
  alpha0(s) = log(q/(1-q));
end
theta0 = [scale0; beta0; alpha0];

objective = @(theta) localNLLGradient(theta, xPref, xProbe, y, nSessions);
options = optimoptions('fminunc', ...
  'Algorithm', 'trust-region', ...
  'SpecifyObjectiveGradient', true, ...
  'Display', char(p.Results.Display), ...
  'MaxFunctionEvaluations', 5e5, ...
  'MaxIterations', 5000, ...
  'OptimalityTolerance', 1e-9, ...
  'FunctionTolerance', 1e-12, ...
  'StepTolerance', 1e-12);

[theta, nll, exitflag, output, gradient, hessian] = fminunc(objective, theta0, options);

theta = full(double(theta));
gradient = full(double(gradient));
hessian = full(double(hessian));
covariance = full(double(pinv(hessian)));
se = sqrt(max(0, diag(covariance)));

scale = theta(1);
betaPref = theta(2:1+nSessions);
alpha = theta(2+nSessions:end);
scaleSE = se(1);

fit = struct();
fit.version = 1;
fit.model = ['P(correct)=0.5+0.5*logistic(alpha_session+' 'betaPref_session*(xPref+scale*xProbe))'];
fit.scale = scale;
fit.scaleSE = scaleSE;
fit.scaleCI95 = scale + 1.96*scaleSE*[-1 1];
fit.betaPrefBySession = betaPref;
fit.betaPrefSEBySession = se(2:1+nSessions);
fit.alphaBySession = alpha;
fit.alphaSEBySession = se(2+nSessions:end);
fit.nSessions = nSessions;
fit.nTrials = sum(nTrialsBySession);
fit.nTrialsBySession = nTrialsBySession;
fit.nCorrectBySession = nCorrectBySession;
fit.sessionFractionCorrect = nCorrectBySession ./ nTrialsBySession;
fit.logLikelihood = -nll;
fit.nParameters = numel(theta);
fit.AIC = 2*fit.nParameters - 2*fit.logLikelihood;
fit.BIC = log(fit.nTrials)*fit.nParameters - 2*fit.logLikelihood;
fit.exitflag = exitflag;
fit.optimizerOutput = output;
fit.gradient = gradient;
fit.gradientInfNorm = norm(gradient, Inf);
fit.hessian = hessian;
fit.covariance = covariance;
fit.fitUsable = isfinite(scale) && isfinite(scaleSE) && ...
  all(isfinite(betaPref)) && all(isfinite(alpha));
end

function [nll, gradient] = localNLLGradient(theta, xPref, xProbe, y, nSessions)
scale = theta(1);
beta = theta(2:1+nSessions);
alpha = theta(2+nSessions:end);

nll = 0;
gScale = 0;
gBeta = zeros(nSessions,1);
gAlpha = zeros(nSessions,1);

for s = 1:nSessions
  z = xPref{s} + scale*xProbe{s};
  eta = alpha(s) + beta(s)*z;
  [sig, p] = logisticAndProbability(eta);
  pSafe = min(max(p, realmin), 1-eps);
  yy = y{s};
  nll = nll - sum(yy.*log(pSafe) + (1-yy).*log(1-pSafe));

  dp = 0.5*sig.*(1-sig);
  dEta = (p-yy).*dp ./ (pSafe.*(1-pSafe));
  gScale = gScale + sum(dEta .* (beta(s)*xProbe{s}));
  gBeta(s) = sum(dEta .* z);
  gAlpha(s) = sum(dEta);
end

gradient = [gScale; gBeta; gAlpha];
end

function [s, p] = logisticAndProbability(eta)
s = zeros(size(eta));
pos = eta >= 0;
s(pos) = 1 ./ (1 + exp(-eta(pos)));
e = exp(eta(~pos));
s(~pos) = e ./ (1 + e);
p = 0.5 + 0.5*s;
end
