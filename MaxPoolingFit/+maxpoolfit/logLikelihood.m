function [ll,pChosen] = logLikelihood(cache, trialIdx, adaptation, noiseSD, cfg, varargin)
%LOGLIKELIHOOD Binary trial-choice log likelihood.
pRDK1 = maxpoolfit.choiceProbability(cache,trialIdx,adaptation,noiseSD,cfg,varargin{:});
choice1 = double(cache.trials.choiceRDK(trialIdx)==1);
pChosen = choice1.*pRDK1 + (1-choice1).*(1-pRDK1);
ll = sum(log(pChosen));
end

