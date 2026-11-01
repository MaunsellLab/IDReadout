function result = profileAdaptationGain(cache,cfg,varargin)
%PROFILEADAPTATIONGAIN Profile prestep-dependent exponential response gain.
%
% The shared subtractive adaptation coefficient may be held fixed while b
% is profiled, or supplied as a second grid dimension, in
%   (Ppre-a*Apre)+exp(-b*Ppre/P0).*(Rstep-Ppre).
% One internal-noise SD is reoptimized for each domain. b=0 is exactly the
% current linear prestep-adaptation model. The adapted prestep baseline is
% preserved; only evidence arriving during the step is gain attenuated.

if ~isfield(cfg,'adaptationExponent'), cfg.adaptationExponent=0; end
if ~isfield(cfg,'adaptationGainStrength'), cfg.adaptationGainStrength=0; end
if ~isfield(cfg,'adaptationReferencePC'), cfg.adaptationReferencePC=40; end
p=inputParser;
p.addParameter('TrialMask',[],@(x)islogical(x)||isnumeric(x));
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('GainGrid',0:0.1:1, ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x))&&all(x>=0));
p.addParameter('Adaptation',cfg.adaptation, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=1);
p.addParameter('AdaptationGrid',[], ...
  @(x)isempty(x)||(isnumeric(x)&&isvector(x)&&all(isfinite(x))&& ...
  all(x>=0)&&all(x<=1)));
p.addParameter('ReferencePC',cfg.adaptationReferencePC, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
assert(string(cache.schemaVersion)==string(cfg.cacheSchemaVersion), ...
  'Cache schema does not match the model configuration.');
if isfield(cache,'kernelStatus')
  assert(string(cache.kernelStatus)=="current", ...
    'Cached temporal predictors are stale; recompute kernels.');
end
if isempty(p.Results.TrialMask)
  mask=cfg.primarySelector(cache.trials);
else
  mask=false(height(cache.trials),1);
  if islogical(p.Results.TrialMask)
    assert(numel(p.Results.TrialMask)==height(cache.trials), ...
      'Logical TrialMask must have one element per cached trial.');
    mask=p.Results.TrialMask(:);
  else
    mask(p.Results.TrialMask)=true;
  end
end
gainGrid=unique(double(p.Results.GainGrid(:)),'stable');
referenceAdaptation=double(p.Results.Adaptation);
if isempty(p.Results.AdaptationGrid)
  adaptationGrid=referenceAdaptation;
else
  adaptationGrid=unique(double(p.Results.AdaptationGrid(:)),'stable');
end
referencePC=double(p.Results.ReferencePC);
domains=unique(string(cache.trials.domain(mask)),'stable');
assert(~isempty(domains),'No trials selected.');

% This profile tests gain compression alone, not its combination with the
% rejected power-weighted subtraction.
cfg=maxpoolfit.configureAdaptationExponent(cfg,0,referencePC);
z1=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,0);
z2=maxpoolfit.noiseDraws(cfg,p.Results.NDraws,1000003);
precomputed=cell(numel(domains),1);
for id=1:numel(domains)
  idx=find(mask & string(cache.trials.domain)==domains(id));
  [B,A,Rpre]=maxpoolfit.deterministicResponse(cache,idx,cfg);
  precomputed{id}=struct('B',B,'A',A,'Rpre',Rpre, ...
    'choice1',double(cache.trials.choiceRDK(idx)==1));
end

nA=numel(adaptationGrid); nB=numel(gainGrid); nD=numel(domains);
nCells=nA*nB;
adaptationColumn=nan(nCells,1); gainColumn=nan(nCells,1);
noiseSDPC=nan(nCells,nD); domainLogLikelihood=nan(nCells,nD);
evalCount=zeros(nCells,nD); row=0;
for ia=1:nA
  adaptation=adaptationGrid(ia);
  for ib=1:nB
    row=row+1; b=gainGrid(ib);
    adaptationColumn(row)=adaptation; gainColumn(row)=b;
    cfgB=maxpoolfit.configureAdaptationGain(cfg,b,referencePC);
    if cfg.verbose
      fprintf('Adaptation %.3g, gain %.3g (%d/%d)\n', ...
        adaptation,b,row,nCells);
    end
    for id=1:nD
      P=precomputed{id};
      objective=@(logSigma)nllAt(exp(logSigma),P,cfgB,adaptation);
      [x,fval,~,out]=fminbnd(objective,log(cfg.noiseBoundsPC(1)), ...
        log(cfg.noiseBoundsPC(2)),optimset('TolX',2e-3,'Display','off'));
      noiseSDPC(row,id)=exp(x);
      domainLogLikelihood(row,id)=-fval;
      evalCount(row,id)=out.funcCount;
    end
  end
end
logLikelihood=sum(domainLogLikelihood,2);
[bestLogLikelihood,bestRow]=max(logLikelihood);
linearRow=find(abs(gainColumn)<=1e-12 & ...
  abs(adaptationColumn-referenceAdaptation)<=1e-12,1);
if isempty(linearRow), deltaLogLFromLinear=nan(nCells,1);
else, deltaLogLFromLinear=logLikelihood-logLikelihood(linearRow); end
profile=table(adaptationColumn,gainColumn,logLikelihood, ...
  logLikelihood-bestLogLikelihood,deltaLogLFromLinear,noiseSDPC, ...
  domainLogLikelihood,evalCount, ...
  'VariableNames',{'adaptation','adaptationGainStrength','logLikelihood', ...
  'deltaLogLFromBest','deltaLogLFromLinear','noiseSDPC', ...
  'domainLogLikelihood','evalCount'});
bestCfg=maxpoolfit.configureAdaptationGain(cfg,gainColumn(bestRow),referencePC);
bestCfg.adaptation=adaptationColumn(bestRow);
linearCfg=maxpoolfit.configureAdaptationGain(cfg,0,referencePC);
result=struct('modelVersion',bestCfg.modelVersion, ...
  'baseModelVersion',linearCfg.modelVersion, ...
  'cacheSchemaVersion',cache.schemaVersion,'domains',domains, ...
  'profile',profile,'bestRow',bestRow, ...
  'adaptationGainStrength',gainColumn(bestRow), ...
  'adaptation',adaptationColumn(bestRow), ...
  'referenceAdaptation',referenceAdaptation, ...
  'adaptationReferencePC',referencePC, ...
  'logLikelihood',bestLogLikelihood, ...
  'noiseSDPC',noiseSDPC(bestRow,:).','nTrials',sum(mask), ...
  'trialMask',mask,'nNoiseDraws',p.Results.NDraws, ...
  'config',bestCfg,'createdAt',datetime('now'));

  function nll=nllAt(sigma,P,cfgB,adaptation)
    p1=maxpoolfit.choiceProbabilityFromResponses( ...
      P.B,P.A,P.Rpre,adaptation,sigma,cfgB,z1,z2,false);
    pChosen=P.choice1.*p1+(1-P.choice1).*(1-p1);
    nll=-sum(log(pChosen));
  end
end
