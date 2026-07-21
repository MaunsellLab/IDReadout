function F = analyzeIDRSignedGains( ...
  D, betaWeibull, lapse, targetPerformance, gainBounds, varargin)
% analyzeIDRSignedGains  Fit signed IDR gains separately by probe offset.
%
% F = analyzeIDRSignedGains(D,beta,lapse,target,gainBounds,...)
%
% Models at each offset:
%   linear            gCP, gCQ, gNP, gNQ
%   noChangeSigned    gCP, gCQ, gNPpos, gNPneg, gNQpos, gNQneg
%   allSigned         positive/negative gains for all four streams
%   opponentNoChange  linear change streams plus opponent-consistent
%                     no-change transforms with a fixed preferred:null ratio
%
% Signed negative predictors retain their physical sign:
%   xPos = max(x,0), xNeg = min(x,0).
% Therefore a positive gNeg makes negative noise reduce effective evidence.
%
% Name-value arguments:
%   PreferredNullRatio  Independently measured ratio (default 2.42)
%   PlotPath            Optional PDF path for signed no-change gains


p = inputParser;
addParameter(p,'PreferredNullRatio',2.42, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p,'PlotPath','',@(x) ischar(x) || isstring(x));
parse(p,varargin{:});
opts = p.Results;

validateInput(D,betaWeibull,lapse,targetPerformance,gainBounds);

T = D.trialTable;
stepFrames = D.stepFrames(:);
XlinearAll = [ ...
  stepMean(D.changePrefNoiseByFrameTrial,stepFrames), ...
  stepMean(D.changeProbeEffectiveNoiseByFrameTrial,stepFrames), ...
  stepMean(D.noChangePrefNoiseByFrameTrial,stepFrames), ...
  stepMean(D.noChangeProbeEffectiveNoiseByFrameTrial,stepFrames)];

offsets = unique(double(T.probeDirDeg),'sorted');
offsets = offsets(isfinite(offsets));
nOffsets = numel(offsets);
emptyResult = struct( ...
  'probeDirDeg',[],'nTrials',[],'nSessions',[], ...
  'linear',[],'noChangeSigned',[],'allSigned',[], ...
  'opponentNoChange',[],'comparisons',[]);
byOffset = repmat(emptyResult,nOffsets,1);

for k = 1:nOffsets
  theta = offsets(k);
  use = double(T.probeDirDeg) == theta;
  X = XlinearAll(use,:);

  common = struct();
  common.stepCoh = double(T.stepDeltaCohPC(use));
  common.correct = double(T.correct(use));
  common.sessionThreshold = double(T.sessionThreshold(use));
  common.betaWeibull = betaWeibull;
  common.lapse = lapse;
  common.targetPerformance = targetPerformance;

  linearStarts = [ ...
     0.95, 1.00, 0.00; ...
     0.40, 0.00, 0.00; ...
    -0.18, 0.00, 0.00; ...
    -0.20, 0.00, 0.00];
  linear = fitGainModel(sprintf('offset_%g_linear',theta), ...
    ["gCP","gCQ","gNP","gNQ"],X,linearStarts,common, ...
    repmat(gainBounds(1),4,1),repmat(gainBounds(2),4,1));

  XncSigned = [X(:,1),X(:,2),posPart(X(:,3)),negPart(X(:,3)), ...
    posPart(X(:,4)),negPart(X(:,4))];
  ncNestedStart = [linear.gain(1);linear.gain(2); ...
    linear.gain(3);linear.gain(3);linear.gain(4);linear.gain(4)];
  ncOpponentLikeStart = [linear.gain(1);linear.gain(2); ...
    -0.3;0.3/opts.PreferredNullRatio;-0.3;0.3/opts.PreferredNullRatio];
  ncStarts = [ncNestedStart,ncOpponentLikeStart,zeros(6,1)];
  noChangeSigned = fitGainModel( ...
    sprintf('offset_%g_noChangeSigned',theta), ...
    ["gCP","gCQ","gNPpos","gNPneg","gNQpos","gNQneg"], ...
    XncSigned,ncStarts,common,repmat(gainBounds(1),6,1), ...
    repmat(gainBounds(2),6,1));

  XallSigned = [posPart(X(:,1)),negPart(X(:,1)), ...
    posPart(X(:,2)),negPart(X(:,2)), ...
    posPart(X(:,3)),negPart(X(:,3)), ...
    posPart(X(:,4)),negPart(X(:,4))];
  allNestedStart = [ ...
    noChangeSigned.gain(1);noChangeSigned.gain(1); ...
    noChangeSigned.gain(2);noChangeSigned.gain(2); ...
    noChangeSigned.gain(3:6)];
  allStarts = [allNestedStart,zeros(8,1)];
  allSigned = fitGainModel(sprintf('offset_%g_allSigned',theta), ...
    ["gCPpos","gCPneg","gCQpos","gCQneg", ...
     "gNPpos","gNPneg","gNQpos","gNQneg"], ...
    XallSigned,allStarts,common,repmat(gainBounds(1),8,1), ...
    repmat(gainBounds(2),8,1));

  r = opts.PreferredNullRatio;
  zNP = -posPart(X(:,3)) + negPart(X(:,3))/r;
  zNQ = -posPart(X(:,4)) + negPart(X(:,4))/r;
  Xopponent = [X(:,1),X(:,2),zNP,zNQ];
  kNP0 = max(0,-linear.gain(3));
  kNQ0 = max(0,-linear.gain(4));
  opponentStarts = [ ...
    linear.gain(1), 1.00, 0.00; ...
    linear.gain(2), 0.30, 0.00; ...
    kNP0,          0.30, 0.00; ...
    kNQ0,          0.30, 0.00];
  opponentLB = [gainBounds(1);gainBounds(1);0;0];
  opponentUB = [gainBounds(2);gainBounds(2); ...
    max(0,gainBounds(2));max(0,gainBounds(2))];
  opponentNoChange = fitGainModel( ...
    sprintf('offset_%g_opponentNoChange',theta), ...
    ["gCP","gCQ","kNPpref","kNQpref"], ...
    Xopponent,opponentStarts,common,opponentLB,opponentUB);

  C = struct();
  C.deltaNLLNoChangeSignedVsLinear = ...
    linear.negLogLikelihood-noChangeSigned.negLogLikelihood;
  C.pNoChangeSignedVsLinear = chiSquareTailFromDeltaNLL( ...
    C.deltaNLLNoChangeSignedVsLinear,2);
  C.deltaNLLAllSignedVsNoChangeSigned = ...
    noChangeSigned.negLogLikelihood-allSigned.negLogLikelihood;
  C.pAllSignedVsNoChangeSigned = chiSquareTailFromDeltaNLL( ...
    C.deltaNLLAllSignedVsNoChangeSigned,2);
  C.deltaNLLOpponentVsFreeNoChangeSigned = ...
    opponentNoChange.negLogLikelihood-noChangeSigned.negLogLikelihood;
  C.pOpponentVsFreeNoChangeSigned = chiSquareTailFromDeltaNLL( ...
    C.deltaNLLOpponentVsFreeNoChangeSigned,2);
  C.deltaAICOpponentMinusFreeNoChangeSigned = ...
    2*(opponentNoChange.negLogLikelihood-noChangeSigned.negLogLikelihood) + ...
    2*(opponentNoChange.nParameters-noChangeSigned.nParameters);

  byOffset(k).probeDirDeg = theta;
  byOffset(k).nTrials = linear.nTrials;
  byOffset(k).nSessions = numel(unique(T.sessionIndex(use)));
  byOffset(k).linear = linear;
  byOffset(k).noChangeSigned = noChangeSigned;
  byOffset(k).allSigned = allSigned;
  byOffset(k).opponentNoChange = opponentNoChange;
  byOffset(k).comparisons = C;
end

F = struct();
F.createdAt = datetime('now');
F.createdBy = mfilename;
F.predictorType = 'rectStep';
F.predictorDefinition = 'signed mean physical coherence over step';
F.probeScaling = 'effective sum over yoked probe streams';
F.preferredNullRatio = opts.PreferredNullRatio;
F.betaWeibull = betaWeibull;
F.lapse = lapse;
F.targetPerformance = targetPerformance;
F.gainBounds = gainBounds(:)';
F.byOffset = byOffset;
F.summaryTable = makeSummaryTable(byOffset,opts.PreferredNullRatio);

if strlength(string(opts.PlotPath)) > 0
  plotSignedNoChangeGains(F,char(opts.PlotPath));
  F.plotPath = string(opts.PlotPath);
end
end

%% ------------------------------------------------------------------------
function validateInput(D,betaWeibull,lapse,targetPerformance,gainBounds)

required = {'trialTable','stepFrames', ...
  'changePrefNoiseByFrameTrial','changeProbeEffectiveNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial','noChangeProbeEffectiveNoiseByFrameTrial'};
for k = 1:numel(required)
  if ~isfield(D,required{k})
    error('analyzeIDRSignedGains:MissingField', ...
      'Input sideGainData is missing %s.',required{k});
  end
end
requiredTableVars = {'probeDirDeg','stepDeltaCohPC','correct', ...
  'sessionThreshold','sessionIndex'};
missing = setdiff(requiredTableVars,D.trialTable.Properties.VariableNames);
if ~isempty(missing)
  error('analyzeIDRSignedGains:MissingTableVariable', ...
    'Input trial table is missing: %s.',strjoin(missing,', '));
end
if ~(isscalar(betaWeibull) && isfinite(betaWeibull) && betaWeibull > 0)
  error('analyzeIDRSignedGains:BadBeta','betaWeibull must be positive.');
end
if ~(isscalar(lapse) && isfinite(lapse) && lapse >= 0 && lapse < 0.5)
  error('analyzeIDRSignedGains:BadLapse','lapse must be in [0,0.5).');
end
if ~(isscalar(targetPerformance) && targetPerformance > 0.5 && targetPerformance < 1)
  error('analyzeIDRSignedGains:BadTarget','targetPerformance must be in (0.5,1).');
end
if ~(isnumeric(gainBounds) && numel(gainBounds) == 2 && ...
    all(isfinite(gainBounds)) && gainBounds(1) < gainBounds(2))
  error('analyzeIDRSignedGains:BadBounds','gainBounds must be [lower upper].');
end
end

%% ------------------------------------------------------------------------
function x = stepMean(noiseByFrameTrial,stepFrames)
x = mean(double(noiseByFrameTrial(stepFrames,:)),1,'omitnan')';
end

function x = posPart(x)
x = max(double(x),0);
end

function x = negPart(x)
x = min(double(x),0);
end

%% ------------------------------------------------------------------------
function fit = fitGainModel(modelName,parameterNames,X,starts,common,lb,ub)

X = double(X);
valid = all(isfinite(X),2) & isfinite(common.stepCoh) & ...
  isfinite(common.correct) & isfinite(common.sessionThreshold);
X = X(valid,:);
stepCoh = common.stepCoh(valid);
correct = common.correct(valid);
sessionThreshold = common.sessionThreshold(valid);
nParameters = size(X,2);
lb = double(lb(:)); ub = double(ub(:));
if numel(lb) ~= nParameters || numel(ub) ~= nParameters || ...
    size(starts,1) ~= nParameters
  error('analyzeIDRSignedGains:ParameterSize', ...
    'Bounds or starts have the wrong size for %s.',modelName);
end

objective = @(g) gainNLL(g,X,stepCoh,correct,sessionThreshold, ...
  common.betaWeibull,common.lapse,common.targetPerformance);
options = optimoptions('fmincon','Display','off','Algorithm','interior-point', ...
  'OptimalityTolerance',1e-8,'StepTolerance',1e-10, ...
  'MaxFunctionEvaluations',15000,'MaxIterations',4000);

bestGain = nan(nParameters,1); bestNLL = inf; bestExitflag = nan;
startNLL = nan(size(starts,2),1); exitflags = nan(size(starts,2),1);
for iStart = 1:size(starts,2)
  x0 = min(max(starts(:,iStart),lb),ub);
  [gain,nll,exitflag] = fmincon( ...
    objective,x0,[],[],[],[],lb,ub,[],options);
  startNLL(iStart)=nll; exitflags(iStart)=exitflag;
  if isfinite(nll) && nll < bestNLL
    bestGain=gain(:); bestNLL=nll; bestExitflag=exitflag;
  end
end
if ~isfinite(bestNLL)
  error('analyzeIDRSignedGains:FitFailure', ...
    'All optimization starts failed for %s.',modelName);
end

H = finiteDifferenceHessian(objective,bestGain);
covariance=nan(nParameters); correlation=nan(nParameters);
SE=nan(nParameters,1); CI95=nan(nParameters,2);
if all(isfinite(H),'all')
  H=(H+H')/2;
  if all(eig(H)>0)
    covariance=inv(H); covariance=(covariance+covariance')/2;
    variance=diag(covariance);
    ok=isfinite(variance) & variance>0;
    SE(ok)=sqrt(variance(ok));
    CI95(ok,:)=bestGain(ok)+[-1 1].*(1.96*SE(ok));
    d=sqrt(diag(covariance));
    correlation=covariance./(d*d');
  end
end

fit=struct();
fit.model=modelName;
fit.parameterNames=string(parameterNames(:))';
fit.gain=bestGain; fit.SE=SE; fit.CI95=CI95;
fit.hessian=H; fit.covariance=covariance; fit.correlation=correlation;
fit.negLogLikelihood=bestNLL; fit.exitflag=bestExitflag;
fit.startNegLogLikelihoods=startNLL; fit.startExitflags=exitflags;
fit.nParameters=nParameters; fit.nTrials=numel(correct);
fit.nEffectiveCohClipped=sum(stepCoh+X*bestGain<0);
end

%% ------------------------------------------------------------------------
function nll = gainNLL(gain,X,stepCoh,correct,sessionThreshold, ...
  betaWeibull,lapse,targetPerformance)
effectiveCoh=max(stepCoh+X*gain(:),0);
p=fixedShapeWeibullP(effectiveCoh,sessionThreshold, ...
  betaWeibull,lapse,targetPerformance);
p=min(max(p,eps),1-eps);
nll=-sum(correct.*log(p)+(1-correct).*log(1-p));
end

function p = fixedShapeWeibullP(coh,threshold,betaWeibull,lapse,targetPerformance)
ratio=(1-lapse-targetPerformance)/(0.5-lapse);
alpha=threshold./((-log(ratio))^(1/betaWeibull));
coh=max(double(coh),0);
p=1-lapse-(0.5-lapse).*exp(-(coh./alpha).^betaWeibull);
end

function H = finiteDifferenceHessian(objective,x)
x=x(:); n=numel(x); H=nan(n); h=3e-4.*max(1,abs(x)); f0=objective(x);
for i=1:n
  ei=zeros(n,1); ei(i)=h(i);
  fp=objective(x+ei); fm=objective(x-ei);
  H(i,i)=(fp-2*f0+fm)/h(i)^2;
  for j=i+1:n
    ej=zeros(n,1); ej(j)=h(j);
    fpp=objective(x+ei+ej); fpm=objective(x+ei-ej);
    fmp=objective(x-ei+ej); fmm=objective(x-ei-ej);
    H(i,j)=(fpp-fpm-fmp+fmm)/(4*h(i)*h(j)); H(j,i)=H(i,j);
  end
end
end

function p = chiSquareTailFromDeltaNLL(deltaNLL,df)
if ~isfinite(deltaNLL) || deltaNLL < -1e-6
  p=nan;
else
  p=gammainc(max(deltaNLL,0),df/2,'upper');
end
end

%% ------------------------------------------------------------------------
function S = makeSummaryTable(byOffset,preferredNullRatio)

n=numel(byOffset);
probeDirDeg=nan(n,1); nTrials=zeros(n,1); nSessions=zeros(n,1);
linearNLL=nan(n,1); noChangeSignedNLL=nan(n,1); allSignedNLL=nan(n,1);
opponentNoChangeNLL=nan(n,1);
deltaNLLNoChangeSigned=nan(n,1); pNoChangeSigned=nan(n,1);
deltaNLLChangeSigned=nan(n,1); pChangeSigned=nan(n,1);
deltaNLLOpponentVsFree=nan(n,1); pOpponentVsFree=nan(n,1);
deltaAICOpponentMinusFree=nan(n,1);
gCP=nan(n,1); seCP=nan(n,1); gCQ=nan(n,1); seCQ=nan(n,1);
gNPpos=nan(n,1); seNPpos=nan(n,1); gNPneg=nan(n,1); seNPneg=nan(n,1);
gNQpos=nan(n,1); seNQpos=nan(n,1); gNQneg=nan(n,1); seNQneg=nan(n,1);
kNPpref=nan(n,1); seKNPpref=nan(n,1); kNQpref=nan(n,1); seKNQpref=nan(n,1);
freeRatioNP=nan(n,1); freeRatioNQ=nan(n,1);
allGCPpos=nan(n,1); allGCPneg=nan(n,1); allGCQpos=nan(n,1); allGCQneg=nan(n,1);

for k=1:n
  R=byOffset(k); N=R.noChangeSigned; A=R.allSigned; O=R.opponentNoChange;
  probeDirDeg(k)=R.probeDirDeg; nTrials(k)=R.nTrials; nSessions(k)=R.nSessions;
  linearNLL(k)=R.linear.negLogLikelihood;
  noChangeSignedNLL(k)=N.negLogLikelihood; allSignedNLL(k)=A.negLogLikelihood;
  opponentNoChangeNLL(k)=O.negLogLikelihood;
  deltaNLLNoChangeSigned(k)=R.comparisons.deltaNLLNoChangeSignedVsLinear;
  pNoChangeSigned(k)=R.comparisons.pNoChangeSignedVsLinear;
  deltaNLLChangeSigned(k)=R.comparisons.deltaNLLAllSignedVsNoChangeSigned;
  pChangeSigned(k)=R.comparisons.pAllSignedVsNoChangeSigned;
  deltaNLLOpponentVsFree(k)=R.comparisons.deltaNLLOpponentVsFreeNoChangeSigned;
  pOpponentVsFree(k)=R.comparisons.pOpponentVsFreeNoChangeSigned;
  deltaAICOpponentMinusFree(k)=R.comparisons.deltaAICOpponentMinusFreeNoChangeSigned;
  gCP(k)=N.gain(1); seCP(k)=N.SE(1); gCQ(k)=N.gain(2); seCQ(k)=N.SE(2);
  gNPpos(k)=N.gain(3); seNPpos(k)=N.SE(3); gNPneg(k)=N.gain(4); seNPneg(k)=N.SE(4);
  gNQpos(k)=N.gain(5); seNQpos(k)=N.SE(5); gNQneg(k)=N.gain(6); seNQneg(k)=N.SE(6);
  if gNPpos(k) < 0 && gNPneg(k) > 0
    freeRatioNP(k)=-gNPpos(k)/gNPneg(k);
  end
  if gNQpos(k) < 0 && gNQneg(k) > 0
    freeRatioNQ(k)=-gNQpos(k)/gNQneg(k);
  end
  kNPpref(k)=O.gain(3); seKNPpref(k)=O.SE(3);
  kNQpref(k)=O.gain(4); seKNQpref(k)=O.SE(4);
  allGCPpos(k)=A.gain(1); allGCPneg(k)=A.gain(2);
  allGCQpos(k)=A.gain(3); allGCQneg(k)=A.gain(4);
end

S=table(probeDirDeg,nTrials,nSessions,preferredNullRatio*ones(n,1), ...
  linearNLL,noChangeSignedNLL,allSignedNLL,opponentNoChangeNLL, ...
  deltaNLLNoChangeSigned,pNoChangeSigned,deltaNLLChangeSigned,pChangeSigned, ...
  deltaNLLOpponentVsFree,pOpponentVsFree,deltaAICOpponentMinusFree, ...
  gCP,seCP,gCQ,seCQ,gNPpos,seNPpos,gNPneg,seNPneg, ...
  gNQpos,seNQpos,gNQneg,seNQneg,freeRatioNP,freeRatioNQ, ...
  kNPpref,seKNPpref,kNQpref,seKNQpref, ...
  allGCPpos,allGCPneg,allGCQpos,allGCQneg, ...
  'VariableNames',{'probeDirDeg','nTrials','nSessions','preferredNullRatio', ...
  'linearNLL','noChangeSignedNLL','allSignedNLL','opponentNoChangeNLL', ...
  'deltaNLLNoChangeSigned','pNoChangeSigned','deltaNLLChangeSigned','pChangeSigned', ...
  'deltaNLLOpponentVsFree','pOpponentVsFree','deltaAICOpponentMinusFree', ...
  'gCP','seCP','gCQ','seCQ','gNPpos','seNPpos','gNPneg','seNPneg', ...
  'gNQpos','seNQpos','gNQneg','seNQneg','freeRatioNP','freeRatioNQ', ...
  'kNPpref','seKNPpref','kNQpref','seKNQpref', ...
  'allGCPpos','allGCPneg','allGCQpos','allGCQneg'});
end

%% ------------------------------------------------------------------------
function plotSignedNoChangeGains(F,pdfPath)

S=F.summaryTable;
fig=figure('Color','w','Visible','off','Units','inches','Position',[1 1 10 7.5]);
tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('IDR signed no-change gains; preferred:null = %.2f', ...
  F.preferredNullRatio),'FontWeight','bold');

ax=nexttile(tl); hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gNPpos,S.seNPpos,'o-','LineWidth',1.2, ...
  'DisplayName','Preferred stream, positive noise');
errorbar(ax,S.probeDirDeg,S.gNPneg,S.seNPneg,'o-','LineWidth',1.2, ...
  'DisplayName','Preferred stream, negative noise');
yline(ax,0,'k:'); ylabel(ax,'Signed gain'); title(ax,'No-change preferred stream');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');

ax=nexttile(tl); hold(ax,'on');
errorbar(ax,S.probeDirDeg,S.gNQpos,S.seNQpos,'o-','LineWidth',1.2, ...
  'DisplayName','Probe stream, positive noise');
errorbar(ax,S.probeDirDeg,S.gNQneg,S.seNQneg,'o-','LineWidth',1.2, ...
  'DisplayName','Probe stream, negative noise');
yline(ax,0,'k:'); xlabel(ax,'Probe offset (deg)'); ylabel(ax,'Signed gain');
title(ax,'No-change probe stream');
legend(ax,'Location','best'); grid(ax,'on'); box(ax,'off');

exportgraphics(fig,pdfPath,'ContentType','vector'); close(fig);
end
