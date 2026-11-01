function result=profileIndependentNoiseMixture(cache,baseFit,varargin)
%PROFILEINDEPENDENTNOISEMIXTURE Sweep independent and response-scaled noise.

p=inputParser;
p.addParameter('IndependentFractionGrid',[0 .05 .1 .2 .4 .7 1], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&& ...
  all(x>=0&x<=1));
p.addParameter('VarianceExponentGrid',0,@(x)isnumeric(x)&&isvector(x)&& ...
  ~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('NDraws',64,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('MinimumGainSessions',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ComputeGainDiagnostics',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('ComputeWinnerDiagnostics',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
required={'trialGroup','trialMask','adaptationByGroup','config'};
assert(all(isfield(baseFit,required)),'A completed design-group fit is required.');
lambda=unique(double(p.Results.IndependentFractionGrid(:)),'stable');
exponent=unique(double(p.Results.VarianceExponentGrid(:)),'stable');
trialGroup=string(baseFit.trialGroup(:));
workCache=cache; workCache.trials.domain=trialGroup;
overall=cell(0,10); groupProfile=table(); gainProfile=table(); winnerProfile=table();
fits=cell(numel(lambda)*numel(exponent),1); fitIndex=0;
for j=1:numel(exponent)
  for k=1:numel(lambda)
    fitIndex=fitIndex+1;
    cfg=maxpoolfit.configureNoiseMixture(baseFit.config,lambda(k),exponent(j));
    cfg.verbose=false;
    fprintf('Independent fraction %g; variance exponent %g (%d/%d)\n', ...
      lambda(k),exponent(j),fitIndex,numel(fits));
    fit=maxpoolfit.fitDesignGroupModel(cache,cfg,trialGroup, ...
      baseFit.adaptationByGroup,'TrialMask',baseFit.trialMask, ...
      'NDraws',p.Results.NDraws);
    fits{fitIndex}=fit;
    G=table(); nConditions=NaN; absoluteRMSE=NaN; ratioRMSE=NaN;
    needPrediction=p.Results.ComputeGainDiagnostics||p.Results.ComputeWinnerDiagnostics;
    if needPrediction
      prediction=maxpoolfit.predictDesignGroupModel(cache,fit,cfg, ...
        'TrialMask',baseFit.trialMask,'NDraws',p.Results.NDraws, ...
        'ReturnWinners',p.Results.ComputeWinnerDiagnostics);
    end
    if p.Results.ComputeGainDiagnostics
      comparison=maxpoolfit.gainComparison(workCache,prediction,prediction, ...
        'TrialMask',baseFit.trialMask);
      G=comparison.summary;
      G=G(ismember(string(G.source),["observed" "primary domain noise"]),:);
      G.independentNoiseFraction=repmat(lambda(k),height(G),1);
      G.noiseVarianceExponent=repmat(exponent(j),height(G),1);
      G.source(string(G.source)=="primary domain noise")="noise-mixture model";
      gainProfile=appendRows(gainProfile,G);
      [nConditions,absoluteRMSE,ratioRMSE]=gainMetrics(G,p.Results.MinimumGainSessions);
    end
    if p.Results.ComputeWinnerDiagnostics
      W=maxpoolfit.winnerFractionDiagnostics(cache,prediction,cfg,trialGroup);
      W.independentNoiseFraction=repmat(lambda(k),height(W),1);
      W.noiseVarianceExponent=repmat(exponent(j),height(W),1);
      winnerProfile=appendRows(winnerProfile,W);
    end
    H=fit.groupFit;
    H.independentNoiseFraction=repmat(lambda(k),height(H),1);
    H.noiseVarianceExponent=repmat(exponent(j),height(H),1);
    groupProfile=appendRows(groupProfile,H);
    overall(end+1,:)={lambda(k),exponent(j),fit.logLikelihood, ... %#ok<AGROW>
      sum(H.evalCount),fit.nTrials,nConditions,absoluteRMSE,ratioRMSE, ...
      mean(H.noiseSDPC),cfg.noiseCorrLengthDeg};
  end
end
overallProfile=cell2table(overall,'VariableNames',{ ...
  'independentNoiseFraction','noiseVarianceExponent','logLikelihood', ...
  'evalCount','nTrials','nGainConditions','absoluteGainRMSE','ratioRMSE', ...
  'meanNoiseSDPC','correlationLengthDeg'});
overallProfile.deltaLogLFromBest=overallProfile.logLikelihood- ...
  max(overallProfile.logLikelihood);
result=struct('overallProfile',overallProfile,'groupProfile',groupProfile, ...
  'gainProfile',gainProfile,'winnerProfile',winnerProfile,'fits',{fits}, ...
  'independentFractionGrid',lambda,'varianceExponentGrid',exponent, ...
  'nNoiseDraws',p.Results.NDraws, ...
  'minimumGainSessions',p.Results.MinimumGainSessions, ...
  'varianceModel',"SD = sigma * (full response / isotropic response)^(p/2); "+ ...
    "variance therefore scales with response^p", ...
  'winnerAssignment',"nearest presented stream direction", ...
  'createdAt',datetime('now'));
end

function [nConditions,absoluteRMSE,ratioRMSE]=gainMetrics(G,minSessions)
O=sortrows(G(string(G.source)=="observed",:), ...
  {'domain','animal','probeOffsetDeg'});
M=sortrows(G(string(G.source)=="noise-mixture model",:), ...
  {'domain','animal','probeOffsetDeg'});
assert(height(O)==height(M),'Observed and predicted gain cells do not align.');
q=O.nSessions>=minSessions; nConditions=sum(q);
if ~any(q), absoluteRMSE=NaN; ratioRMSE=NaN; return; end
d=M.changedDriftGain(q)-O.changedDriftGain(q);
r=M.changedProbeGain(q)-O.changedProbeGain(q);
absoluteRMSE=sqrt(mean([d.^2;r.^2]));
ratioRMSE=sqrt(mean((M.changedProbeOverDrift(q)- ...
  O.changedProbeOverDrift(q)).^2));
end

function A=appendRows(A,B)
if width(A)==0, A=B; else, A=[A;B]; end
end
