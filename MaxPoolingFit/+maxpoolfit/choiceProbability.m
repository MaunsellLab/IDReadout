function [pRDK1, details] = choiceProbability(cache, trialIdx, adaptation, noiseSD, cfg, varargin)
%CHOICEPROBABILITY Integrate hard-max choices over internal Gaussian noise.
p = inputParser;
p.addParameter('NDraws',cfg.nNoiseDraws,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
p.addParameter('ReturnWinners',false,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
nDraws = p.Results.NDraws;

z1 = maxpoolfit.noiseDraws(cfg,nDraws,0);
z2 = maxpoolfit.noiseDraws(cfg,nDraws,1000003);
[B,A,P] = maxpoolfit.deterministicResponse(cache,trialIdx,cfg);
[pRDK1,details] = maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,P,adaptation,noiseSD,cfg,z1,z2,p.Results.ReturnWinners);
end
