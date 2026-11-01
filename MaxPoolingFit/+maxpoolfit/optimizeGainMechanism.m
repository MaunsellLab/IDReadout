function result=optimizeGainMechanism(cache,baseFit,varargin)
%OPTIMIZEGAINMECHANISM Jointly fit mechanism parameters to changed-side gains.
%
% Effective directional width, prestep adaptation, and directional-noise
% correlation length are optimized together. At every candidate, one
% upstream noise SD per target design group is profiled by trial-level
% choice likelihood. The outer objective gives equal weight to changed-side
% drift and probe gains in every retained animal-design-offset cell.

p=inputParser;
p.addParameter('TargetGroups',["IDQ-randomized" "IDS"], ...
  @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('DirectionalWidthBoundsDeg',[10 37.5],@validBounds);
p.addParameter('AdaptationBounds',[0 .8],@validBounds);
p.addParameter('CorrelationLengthBoundsDeg',[20 180],@validBounds);
p.addParameter('InitialPoint',[20 .4 120], ...
  @(x)isnumeric(x)&&numel(x)==3&&all(isfinite(x)));
p.addParameter('NInitialPoints',12,@positiveInteger);
p.addParameter('NStarts',3,@positiveInteger);
p.addParameter('MaxFunctionEvaluationsPerStart',30,@positiveInteger);
p.addParameter('NDraws',64,@evenDrawCount);
p.addParameter('FinalNDraws',256,@evenDrawCount);
p.addParameter('MinimumGainSessions',3,@positiveInteger);
p.addParameter('SearchSeed',9187,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
p.addParameter('CheckpointPath',"",@(x)isstring(x)||ischar(x));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
required={'trialGroup','trialMask','adaptationByGroup','config','modelVersion'};
assert(all(isfield(baseFit,required)), ...
  'A completed design-group fit is required.');
assert(numel(baseFit.trialGroup)==height(cache.trials), ...
  'Design-group labels do not match the cache.');

if ischar(p.Results.TargetGroups), targetGroups=string(p.Results.TargetGroups);
else, targetGroups=unique(string(p.Results.TargetGroups(:)),'stable'); end
trialGroup=string(baseFit.trialGroup(:));
targetMask=logical(baseFit.trialMask(:)) & ismember(trialGroup,targetGroups);
assert(all(ismember(targetGroups,unique(trialGroup(targetMask),'stable'))), ...
  'One or more TargetGroups have no selected primary trials.');
workCache=cache; workCache.trials.domain=trialGroup;

lb=[p.Results.DirectionalWidthBoundsDeg(1) ...
  p.Results.AdaptationBounds(1) p.Results.CorrelationLengthBoundsDeg(1)];
ub=[p.Results.DirectionalWidthBoundsDeg(2) ...
  p.Results.AdaptationBounds(2) p.Results.CorrelationLengthBoundsDeg(2)];
x0=double(p.Results.InitialPoint(:)');
assert(all(x0>=lb & x0<=ub),'InitialPoint lies outside the parameter bounds.');

evaluation=table(); groupEvaluation=table(); gainEvaluation=table();
evaluationCount=0;
initialPoints=spaceFillingPoints(p.Results.NInitialPoints,lb,ub, ...
  p.Results.SearchSeed);
initialPoints(1,:)=x0;

for k=1:size(initialPoints,1)
  evaluate(initialPoints(k,:),p.Results.NDraws,"initial");
end
[~,order]=sort(evaluation.changedGainRMSE);
nStarts=min(p.Results.NStarts,height(evaluation));
startPoints=evaluation{order(1:nStarts), ...
  {'directionalWidthDeg','adaptation','correlationLengthDeg'}};
options=optimset('Display','off','TolX',2e-2,'TolFun',2e-3, ...
  'MaxFunEvals',p.Results.MaxFunctionEvaluationsPerStart, ...
  'MaxIter',p.Results.MaxFunctionEvaluationsPerStart);
for k=1:nStarts
  fprintf('Joint gain refinement %d/%d\n',k,nStarts);
  y0=toUnbounded(startPoints(k,:),lb,ub);
  fminsearch(@refinementObjective,y0,options);
end

[~,gainRow]=min(evaluation.changedGainRMSE);
[~,choiceRow]=max(evaluation.targetLogLikelihood);
evaluation.isPareto=paretoRows(evaluation.changedGainRMSE, ...
  -evaluation.targetLogLikelihood);
gainPoint=evaluation{gainRow, ...
  {'directionalWidthDeg','adaptation','correlationLengthDeg'}};
choicePoint=evaluation{choiceRow, ...
  {'directionalWidthDeg','adaptation','correlationLengthDeg'}};

gainFinal=evaluateFinal(gainPoint,p.Results.FinalNDraws);
if isequal(gainPoint,choicePoint), choiceFinal=gainFinal;
else, choiceFinal=evaluateFinal(choicePoint,p.Results.FinalNDraws); end

result=struct('evaluation',evaluation,'groupEvaluation',groupEvaluation, ...
  'gainEvaluation',gainEvaluation,'paretoEvaluation',evaluation(evaluation.isPareto,:), ...
  'initialPoints',initialPoints,'targetGroups',targetGroups, ...
  'trialMask',targetMask,'parameterLowerBounds',lb,'parameterUpperBounds',ub, ...
  'nNoiseDraws',p.Results.NDraws,'finalNNoiseDraws',p.Results.FinalNDraws, ...
  'minimumGainSessions',p.Results.MinimumGainSessions, ...
  'gainOptimalPoint',gainPoint,'choiceOptimalPoint',choicePoint, ...
  'gainOptimalFit',gainFinal.fit,'gainOptimalPrediction',gainFinal.prediction, ...
  'gainOptimalComparison',gainFinal.comparison, ...
  'gainOptimalMetrics',gainFinal.metrics, ...
  'choiceOptimalFit',choiceFinal.fit, ...
  'choiceOptimalPrediction',choiceFinal.prediction, ...
  'choiceOptimalComparison',choiceFinal.comparison, ...
  'choiceOptimalMetrics',choiceFinal.metrics, ...
  'gainObjective',"equal weight per retained animal-design-offset cell; "+ ...
    "equal weight for changed drift and changed probe gains", ...
  'choiceNuisanceFit',"one upstream noise SD per target design group", ...
  'optimizer',"space-filling initialization plus multistart bounded fminsearch", ...
  'createdAt',datetime('now'));

  function value=refinementObjective(y)
    row=evaluate(fromUnbounded(y,lb,ub),p.Results.NDraws,"refinement");
    value=row.changedGainRMSE;
  end

  function row=evaluate(x,nDraws,stage)
    x=double(x(:)');
    evaluationCount=evaluationCount+1;
    fprintf(['Evaluation %d: sigma %.3g deg; adaptation %.3g; ' ...
      'correlation %.3g deg\n'],evaluationCount,x(1),x(2),x(3));
    out=fitPoint(x,nDraws);
    M=out.metrics; H=out.fit.groupFit;
    row=table(evaluationCount,string(stage),x(1),x(2),x(3), ...
      out.fit.logLikelihood,sum(H.logLikelihood),sum(H.evalCount), ...
      out.fit.nTrials,M.nConditions,M.coefficientRMSE,M.gainRMSE, ...
      M.driftGainRMSE,M.probeGainRMSE,M.ratioRMSE,mean(H.noiseSDPC), ...
      'VariableNames',{'evaluation','stage','directionalWidthDeg', ...
      'adaptation','correlationLengthDeg','logLikelihood', ...
      'targetLogLikelihood','noiseEvalCount','nTrials','nGainConditions', ...
      'changedCoefficientRMSE','changedGainRMSE','changedDriftGainRMSE', ...
      'changedProbeGainRMSE','ratioRMSE','meanTargetNoiseSDPC'});
    evaluation=appendRows(evaluation,row);
    H.evaluation=repmat(evaluationCount,height(H),1);
    H.directionalWidthDeg=repmat(x(1),height(H),1);
    H.correlationLengthDeg=repmat(x(3),height(H),1);
    groupEvaluation=appendRows(groupEvaluation,H);
    G=out.gains; G.evaluation=repmat(evaluationCount,height(G),1);
    G.directionalWidthDeg=repmat(x(1),height(G),1);
    G.adaptation=repmat(x(2),height(G),1);
    G.correlationLengthDeg=repmat(x(3),height(G),1);
    gainEvaluation=appendRows(gainEvaluation,G);
    saveCheckpoint();
  end

  function out=evaluateFinal(x,nDraws)
    fprintf(['Final evaluation: sigma %.3g deg; adaptation %.3g; ' ...
      'correlation %.3g deg; %d draws\n'],x(1),x(2),x(3),nDraws);
    out=fitPoint(x,nDraws);
  end

  function out=fitPoint(x,nDraws)
    cfg=baseFit.config;
    cfg.sigmaMTDeg=x(1); cfg.sigmaAdaptDeg=x(1);
    cfg.noiseCorrLengthDeg=x(3); cfg.verbose=false;
    cfg.modelVersion=string(baseFit.modelVersion)+"-joint-gain-mechanism-fit";
    adaptationByGroup=table(targetGroups,repmat(x(2),numel(targetGroups),1), ...
      'VariableNames',{'analysisGroup','adaptation'});
    fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
      adaptationByGroup,'TrialMask',targetMask,'NDraws',nDraws);
    prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg, ...
      'TrialMask',targetMask,'NDraws',nDraws);
    comparison=maxpoolfit.gainComparison(workCache,prediction,prediction, ...
      'TrialMask',targetMask);
    G=comparison.summary;
    G=G(ismember(string(G.source),["observed" "primary domain noise"]),:);
    G.source(string(G.source)=="primary domain noise")="joint gain model";
    metrics=gainMetrics(G,targetGroups,p.Results.MinimumGainSessions);
    out=struct('fit',fit,'prediction',prediction, ...
      'comparison',comparison,'gains',G,'metrics',metrics);
  end

  function saveCheckpoint()
    checkpointPath=string(p.Results.CheckpointPath);
    if strlength(checkpointPath)==0, return; end
    checkpoint=struct('evaluation',evaluation, ...
      'groupEvaluation',groupEvaluation,'gainEvaluation',gainEvaluation, ...
      'targetGroups',targetGroups,'parameterLowerBounds',lb, ...
      'parameterUpperBounds',ub,'nNoiseDraws',p.Results.NDraws, ...
      'minimumGainSessions',p.Results.MinimumGainSessions, ...
      'updatedAt',datetime('now')); %#ok<NASGU>
    folder=fileparts(checkpointPath);
    if strlength(folder)>0 && ~isfolder(folder), mkdir(folder); end
    temporaryPath=checkpointPath+".tmp.mat";
    save(temporaryPath,'checkpoint','-v7.3');
    movefile(temporaryPath,checkpointPath,'f');
  end
end

function M=gainMetrics(G,targetGroups,minSessions)
O=sortrows(G(ismember(string(G.domain),targetGroups)& ...
  string(G.source)=="observed",:),{'domain','animal','probeOffsetDeg'});
P=sortrows(G(ismember(string(G.domain),targetGroups)& ...
  string(G.source)=="joint gain model",:),{'domain','animal','probeOffsetDeg'});
assert(height(O)==height(P)&&isequal(string(O.domain),string(P.domain))&& ...
  isequal(string(O.animal),string(P.animal))&& ...
  isequal(O.probeOffsetDeg,P.probeOffsetDeg), ...
  'Observed and predicted target gain cells do not align.');
q=O.nSessions>=minSessions;
assert(any(q),'No target gain cells meet MinimumGainSessions.');
driftGain=P.changedDriftGain(q)-O.changedDriftGain(q);
probeGain=P.changedProbeGain(q)-O.changedProbeGain(q);
driftCoefficient=P.changedDriftCoefficientPerPC(q)- ...
  O.changedDriftCoefficientPerPC(q);
probeCoefficient=P.changedProbeCoefficientPerPC(q)- ...
  O.changedProbeCoefficientPerPC(q);
ratio=P.changedProbeOverDrift(q)-O.changedProbeOverDrift(q);
M=struct('nConditions',sum(q), ...
  'coefficientRMSE',sqrt(mean([driftCoefficient.^2;probeCoefficient.^2])), ...
  'gainRMSE',sqrt(mean([driftGain.^2;probeGain.^2])), ...
  'driftGainRMSE',sqrt(mean(driftGain.^2)), ...
  'probeGainRMSE',sqrt(mean(probeGain.^2)), ...
  'ratioRMSE',sqrt(mean(ratio.^2)));
end

function X=spaceFillingPoints(n,lb,ub,seed)
old=rng; cleanup=onCleanup(@()rng(old)); %#ok<NASGU>
rng(seed,'twister');
if n==1, U=.5*ones(1,numel(lb));
elseif exist('lhsdesign','file')==2
  U=lhsdesign(n,numel(lb),'criterion','maximin','iterations',20);
else
  U=zeros(n,numel(lb));
  for j=1:numel(lb), U(:,j)=(randperm(n)'-.5)/n; end
end
X=lb+(ub-lb).*U;
end

function y=toUnbounded(x,lb,ub)
u=(x-lb)./(ub-lb); u=min(max(u,1e-6),1-1e-6);
y=log(u./(1-u));
end

function x=fromUnbounded(y,lb,ub)
u=1./(1+exp(-double(y(:)'))); x=lb+(ub-lb).*u;
end

function q=paretoRows(gainError,nll)
n=numel(gainError); q=true(n,1);
for i=1:n
  q(i)=~any((gainError<=gainError(i) & nll<=nll(i)) & ...
    (gainError<gainError(i) | nll<nll(i)));
end
end

function tf=validBounds(x)
tf=isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)<x(2);
end
function tf=positiveInteger(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&mod(x,1)==0;
end
function tf=evenDrawCount(x)
tf=positiveInteger(x)&&x>=2&&mod(x,2)==0;
end
function A=appendRows(A,B)
if width(A)==0, A=B; else, A=[A;B]; end
end
