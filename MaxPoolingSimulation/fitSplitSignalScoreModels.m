function fit = fitSplitSignalScoreModels(outcome,scores,cfg)
% FITSPLITSIGNALSCOREMODELS  Refit decision scale and max softness.

fit=struct(); fit.logLikelihood=-inf(numel(cfg.PoolingModels),1);
fit.detail=cell(numel(cfg.PoolingModels),1); fit.bestTemperature=nan;
for m=1:numel(cfg.PoolingModels)
    name=cfg.PoolingModels{m};
    if strcmpi(name,'max')
        for k=1:size(scores.max,2)
            candidate=fitLogisticIRLS(scores.max(:,k),outcome);
            if candidate.logLikelihood>fit.logLikelihood(m)
                fit.logLikelihood(m)=candidate.logLikelihood;
                fit.detail{m}=candidate;
                fit.bestTemperature=cfg.SoftmaxTemperatures(k);
            end
        end
    else
        candidate=fitLogisticIRLS(scores.(name),outcome);
        fit.logLikelihood(m)=candidate.logLikelihood; fit.detail{m}=candidate;
    end
end
% Every family fits an intercept and decision scale.  Profiling max softness
% adds one parameter, which must be penalized to avoid favoring max merely
% because it spans a larger family.
n=max(numel(outcome),1);
fit.parameterCount=2*ones(numel(cfg.PoolingModels),1);
maxIndex=find(strcmpi(cfg.PoolingModels,'max'),1);
if ~isempty(maxIndex), fit.parameterCount(maxIndex)=3; end
switch lower(cfg.SelectionCriterion)
    case 'bic'
        fit.selectionScore=fit.logLikelihood-0.5*fit.parameterCount*log(n);
    case 'aic'
        fit.selectionScore=fit.logLikelihood-fit.parameterCount;
    case {'loglikelihood','ll'}
        fit.selectionScore=fit.logLikelihood;
    otherwise
        error('Unknown SelectionCriterion: %s',cfg.SelectionCriterion);
end
[~,fit.winnerIndex]=max(fit.selectionScore);
end
