function F = analyzeIDRMatchedGains( ...
  D, betaWeibull, lapse, targetPerformance, gainBounds, varargin)
% analyzeIDRMatchedGains  Fit absolute IDR side-by-stream gains by offset.
%
% F = analyzeIDRMatchedGains(D,beta,lapse,target,gainBounds,...)
%
% For each probe offset theta, jointly fit rectangular step predictors:
%
%   cEff = cStep + gCP*nCP + gCQ*nCQ + gNP*nNP + gNQ*nNQ
%
% where C/N denote change/no-change patch and P/Q denote preferred/probe
% stream. Probe predictors are the effective sum over yoked streams, making
% gCQ and gNQ the common per-stream gain under an equal-gain linear model.
% Every offset is fit independently. Physical coherence signs are retained.
%
% A nested preferred-only model sets gCQ=gNQ=0. The full fit retains its
% complete covariance and correlation matrices in F.byOffset.
%
% Name-value argument:
%   PlotPath   Optional output path for an absolute-gain PDF.


p = inputParser;
addParameter(p, 'PlotPath', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});
opts = p.Results;

validateInput(D, betaWeibull, lapse, targetPerformance, gainBounds);

T = D.trialTable;
stepFrames = D.stepFrames(:);
Xall = [ ...
  stepMean(D.changePrefNoiseByFrameTrial, stepFrames), ...
  stepMean(D.changeProbeEffectiveNoiseByFrameTrial, stepFrames), ...
  stepMean(D.noChangePrefNoiseByFrameTrial, stepFrames), ...
  stepMean(D.noChangeProbeEffectiveNoiseByFrameTrial, stepFrames)];

offsets = unique(double(T.probeDirDeg), 'sorted');
offsets = offsets(isfinite(offsets));
nOffsets = numel(offsets);

emptyResult = struct( ...
  'probeDirDeg', [], 'nTrials', [], 'nSessions', [], ...
  'full', [], 'preferredOnly', [], 'deltaNLLProbeTerms', [], ...
  'lrtPProbeTerms', [], 'deltaAICPreferredMinusFull', []);
byOffset = repmat(emptyResult, nOffsets, 1);

for k = 1:nOffsets
  theta = offsets(k);
  use = double(T.probeDirDeg) == theta;

  common = struct();
  common.stepCoh = double(T.stepDeltaCohPC(use));
  common.correct = double(T.correct(use));
  common.sessionThreshold = double(T.sessionThreshold(use));
  common.betaWeibull = betaWeibull;
  common.lapse = lapse;
  common.targetPerformance = targetPerformance;

  X = Xall(use,:);
  preferredStarts = [0.95, 1.00, 0.00; -0.18, 0.00, 0.00];
  preferredOnly = fitGainModel( ...
    sprintf('offset_%g_preferredOnly',theta), ...
    ["gCP","gNP"], X(:,[1 3]), preferredStarts, common, gainBounds);

  % The first full-model start is the fitted nested model itself. This
  % ensures that the unrestricted optimizer starts at a solution with the
  % preferred-only NLL before the two probe gains are released.
  fullStarts = [ ...
     preferredOnly.gain(1), 0.95, 1.00, 0.00; ...
     0.00,                  0.40, 0.00, 0.00; ...
     preferredOnly.gain(2),-0.18, 0.00, 0.00; ...
     0.00,                 -0.10, 0.00, 0.00];
  full = fitGainModel( ...
    sprintf('offset_%g_fourGain',theta), ...
    ["gCP","gCQ","gNP","gNQ"], X, fullStarts, common, gainBounds);

  deltaNLL = preferredOnly.negLogLikelihood - full.negLogLikelihood;
  % LRT statistic is 2*deltaNLL with 2 df, whose upper-tail probability is
  % exp(-deltaNLL). Clamp tiny negative numerical differences to zero.
  lrtP = exp(-max(deltaNLL,0));
  deltaAIC = 2*deltaNLL - 2*(full.nParameters-preferredOnly.nParameters);

  byOffset(k).probeDirDeg = theta;
  byOffset(k).nTrials = full.nTrials;
  byOffset(k).nSessions = numel(unique(T.sessionIndex(use)));
  byOffset(k).full = full;
  byOffset(k).preferredOnly = preferredOnly;
  byOffset(k).deltaNLLProbeTerms = deltaNLL;
  byOffset(k).lrtPProbeTerms = lrtP;
  byOffset(k).deltaAICPreferredMinusFull = deltaAIC;
end

summaryTable = makeSummaryTable(byOffset);

F = struct();
F.createdAt = datetime('now');
F.createdBy = mfilename;
F.predictorType = 'rectStep';
F.predictorDefinition = 'signed mean physical coherence over step';
F.probeScaling = 'effective sum over yoked probe streams';
F.modelEquation = 'cEff=cStep+gCP*nCP+gCQ*nCQ+gNP*nNP+gNQ*nNQ';
F.parameterNames = ["gCP","gCQ","gNP","gNQ"];
F.betaWeibull = betaWeibull;
F.lapse = lapse;
F.targetPerformance = targetPerformance;
F.gainBounds = gainBounds(:)';
F.byOffset = byOffset;
F.summaryTable = summaryTable;

if strlength(string(opts.PlotPath)) > 0
  plotAbsoluteGains(F, char(opts.PlotPath));
  F.plotPath = string(opts.PlotPath);
end
end

%% ------------------------------------------------------------------------
function validateInput(D, betaWeibull, lapse, targetPerformance, gainBounds)

required = { ...
  'trialTable','stepFrames', ...
  'changePrefNoiseByFrameTrial','changeProbeEffectiveNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial','noChangeProbeEffectiveNoiseByFrameTrial'};
for k = 1:numel(required)
  if ~isfield(D,required{k})
    error('analyzeIDRMatchedGains:MissingField', ...
      'Input sideGainData is missing %s.', required{k});
  end
end

requiredTableVars = { ...
  'probeDirDeg','stepDeltaCohPC','correct','sessionThreshold','sessionIndex'};
missing = setdiff(requiredTableVars,D.trialTable.Properties.VariableNames);
if ~isempty(missing)
  error('analyzeIDRMatchedGains:MissingTableVariable', ...
    'Input trial table is missing: %s.', strjoin(missing,', '));
end

if ~(isscalar(betaWeibull) && isfinite(betaWeibull) && betaWeibull > 0)
  error('analyzeIDRMatchedGains:BadBeta','betaWeibull must be positive.');
end
if ~(isscalar(lapse) && isfinite(lapse) && lapse >= 0 && lapse < 0.5)
  error('analyzeIDRMatchedGains:BadLapse','lapse must be in [0,0.5).');
end
if ~(isscalar(targetPerformance) && targetPerformance > 0.5 && targetPerformance < 1)
  error('analyzeIDRMatchedGains:BadTarget','targetPerformance must be in (0.5,1).');
end
if ~(isnumeric(gainBounds) && numel(gainBounds) == 2 && ...
    all(isfinite(gainBounds)) && gainBounds(1) < gainBounds(2))
  error('analyzeIDRMatchedGains:BadBounds','gainBounds must be [lower upper].');
end
end

%% ------------------------------------------------------------------------
function x = stepMean(noiseByFrameTrial, stepFrames)

x = mean(double(noiseByFrameTrial(stepFrames,:)),1,'omitnan')';
end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName, parameterNames, X, starts, common, gainBounds)

X = double(X);
valid = all(isfinite(X),2) & isfinite(common.stepCoh) & ...
  isfinite(common.correct) & isfinite(common.sessionThreshold);
X = X(valid,:);
stepCoh = common.stepCoh(valid);
correct = common.correct(valid);
sessionThreshold = common.sessionThreshold(valid);

nParameters = size(X,2);
if size(starts,1) ~= nParameters
  error('analyzeIDRMatchedGains:StartSize', ...
    'Start matrix for %s has the wrong number of rows.', modelName);
end

lb = repmat(gainBounds(1),nParameters,1);
ub = repmat(gainBounds(2),nParameters,1);
objective = @(g) gainNLL(g,X,stepCoh,correct,sessionThreshold, ...
  common.betaWeibull,common.lapse,common.targetPerformance);

options = optimoptions('fmincon', ...
  'Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',10000,'MaxIterations',3000);

bestGain = nan(nParameters,1);
bestNLL = inf;
bestExitflag = nan;
startNLL = nan(size(starts,2),1);
exitflags = nan(size(starts,2),1);
for iStart = 1:size(starts,2)
  x0 = min(max(starts(:,iStart),lb),ub);
  [gain,nll,exitflag] = fmincon( ...
    objective,x0,[],[],[],[],lb,ub,[],options);
  startNLL(iStart) = nll;
  exitflags(iStart) = exitflag;
  if isfinite(nll) && nll < bestNLL
    bestGain = gain(:);
    bestNLL = nll;
    bestExitflag = exitflag;
  end
end

if ~isfinite(bestNLL)
  error('analyzeIDRMatchedGains:FitFailure', ...
    'All optimization starts failed for %s.',modelName);
end

H = finiteDifferenceHessian(objective,bestGain);
covariance = nan(nParameters);
correlation = nan(nParameters);
SE = nan(nParameters,1);
CI95 = nan(nParameters,2);
if all(isfinite(H),'all')
  H = (H+H')/2;
  if all(eig(H) > 0)
    covariance = inv(H);
    covariance = (covariance+covariance')/2;
    variance = diag(covariance);
    ok = isfinite(variance) & variance > 0;
    SE(ok) = sqrt(variance(ok));
    CI95(ok,:) = bestGain(ok) + [-1 1].*(1.96*SE(ok));
    d = sqrt(diag(covariance));
    correlation = covariance./(d*d');
  end
end

fit = struct();
fit.model = modelName;
fit.parameterNames = string(parameterNames(:))';
fit.gain = bestGain;
fit.SE = SE;
fit.CI95 = CI95;
fit.hessian = H;
fit.covariance = covariance;
fit.correlation = correlation;
fit.negLogLikelihood = bestNLL;
fit.exitflag = bestExitflag;
fit.startNegLogLikelihoods = startNLL;
fit.startExitflags = exitflags;
fit.nParameters = nParameters;
fit.nTrials = numel(correct);
fit.nEffectiveCohClipped = sum(stepCoh + X*bestGain < 0);
end

%% ------------------------------------------------------------------------
function nll = gainNLL(gain,X,stepCoh,correct,sessionThreshold, ...
  betaWeibull,lapse,targetPerformance)

effectiveCoh = max(stepCoh + X*gain(:),0);
p = fixedShapeWeibullP(effectiveCoh,sessionThreshold, ...
  betaWeibull,lapse,targetPerformance);
p = min(max(p,eps),1-eps);
nll = -sum(correct.*log(p) + (1-correct).*log(1-p));
end

%% ------------------------------------------------------------------------
function p = fixedShapeWeibullP(coh,threshold,betaWeibull,lapse,targetPerformance)

ratio = (1-lapse-targetPerformance)/(0.5-lapse);
alpha = threshold./((-log(ratio))^(1/betaWeibull));
coh = max(double(coh),0);
p = 1-lapse-(0.5-lapse).*exp(-(coh./alpha).^betaWeibull);
end

%% ------------------------------------------------------------------------
function H = finiteDifferenceHessian(objective,x)

x = x(:);
n = numel(x);
H = nan(n);
h = 3e-4.*max(1,abs(x));
f0 = objective(x);
for i = 1:n
  ei = zeros(n,1);
  ei(i) = h(i);
  fp = objective(x+ei);
  fm = objective(x-ei);
  H(i,i) = (fp-2*f0+fm)/h(i)^2;
  for j = i+1:n
    ej = zeros(n,1);
    ej(j) = h(j);
    fpp = objective(x+ei+ej);
    fpm = objective(x+ei-ej);
    fmp = objective(x-ei+ej);
    fmm = objective(x-ei-ej);
    H(i,j) = (fpp-fpm-fmp+fmm)/(4*h(i)*h(j));
    H(j,i) = H(i,j);
  end
end
end

%% ------------------------------------------------------------------------
function S = makeSummaryTable(byOffset)

n = numel(byOffset);
probeDirDeg = nan(n,1);
nTrials = zeros(n,1);
nSessions = zeros(n,1);
fullNLL = nan(n,1);
preferredOnlyNLL = nan(n,1);
deltaNLLProbeTerms = nan(n,1);
lrtPProbeTerms = nan(n,1);
deltaAICPreferredMinusFull = nan(n,1);
gCP = nan(n,1); seCP = nan(n,1); ciCPlo = nan(n,1); ciCPhi = nan(n,1);
gCQ = nan(n,1); seCQ = nan(n,1); ciCQlo = nan(n,1); ciCQhi = nan(n,1);
gNP = nan(n,1); seNP = nan(n,1); ciNPlo = nan(n,1); ciNPhi = nan(n,1);
gNQ = nan(n,1); seNQ = nan(n,1); ciNQlo = nan(n,1); ciNQhi = nan(n,1);
fullExitflag = nan(n,1);
nEffectiveCohClipped = zeros(n,1);

for k = 1:n
  R = byOffset(k);
  U = R.full;
  probeDirDeg(k) = R.probeDirDeg;
  nTrials(k) = R.nTrials;
  nSessions(k) = R.nSessions;
  fullNLL(k) = U.negLogLikelihood;
  preferredOnlyNLL(k) = R.preferredOnly.negLogLikelihood;
  deltaNLLProbeTerms(k) = R.deltaNLLProbeTerms;
  lrtPProbeTerms(k) = R.lrtPProbeTerms;
  deltaAICPreferredMinusFull(k) = R.deltaAICPreferredMinusFull;
  gCP(k)=U.gain(1); seCP(k)=U.SE(1); ciCPlo(k)=U.CI95(1,1); ciCPhi(k)=U.CI95(1,2);
  gCQ(k)=U.gain(2); seCQ(k)=U.SE(2); ciCQlo(k)=U.CI95(2,1); ciCQhi(k)=U.CI95(2,2);
  gNP(k)=U.gain(3); seNP(k)=U.SE(3); ciNPlo(k)=U.CI95(3,1); ciNPhi(k)=U.CI95(3,2);
  gNQ(k)=U.gain(4); seNQ(k)=U.SE(4); ciNQlo(k)=U.CI95(4,1); ciNQhi(k)=U.CI95(4,2);
  fullExitflag(k) = U.exitflag;
  nEffectiveCohClipped(k) = U.nEffectiveCohClipped;
end

S = table(probeDirDeg,nTrials,nSessions,fullNLL,preferredOnlyNLL, ...
  deltaNLLProbeTerms,lrtPProbeTerms,deltaAICPreferredMinusFull, ...
  gCP,seCP,ciCPlo,ciCPhi,gCQ,seCQ,ciCQlo,ciCQhi, ...
  gNP,seNP,ciNPlo,ciNPhi,gNQ,seNQ,ciNQlo,ciNQhi, ...
  fullExitflag,nEffectiveCohClipped);
end

%% ------------------------------------------------------------------------
function plotAbsoluteGains(F,pdfPath)

S = F.summaryTable;
fig = figure('Color','w','Visible','off', ...
  'Units','inches','Position',[1 1 10 7.5]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
title(tl,'IDR absolute rectangular-step gains by probe offset', ...
  'FontWeight','bold');

ax = nexttile(tl);
hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gCP,S.seCP,'o-','LineWidth',1.2, ...
  'DisplayName','Change preferred');
errorbar(ax,S.probeDirDeg,S.gCQ,S.seCQ,'o-','LineWidth',1.2, ...
  'DisplayName','Change probe');
yline(ax,0,'k:');
ylabel(ax,'Absolute gain');
title(ax,'Change patch');
legend(ax,'Location','best');
grid(ax,'on'); box(ax,'off');

ax = nexttile(tl);
hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gNP,S.seNP,'o-','LineWidth',1.2, ...
  'DisplayName','No-change preferred');
errorbar(ax,S.probeDirDeg,S.gNQ,S.seNQ,'o-','LineWidth',1.2, ...
  'DisplayName','No-change probe');
yline(ax,0,'k:');
xlabel(ax,'Probe offset (deg)');
ylabel(ax,'Absolute gain');
title(ax,'No-change patch');
legend(ax,'Location','best');
grid(ax,'on'); box(ax,'off');

linkaxes(findall(fig,'Type','axes'),'x');
exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end
