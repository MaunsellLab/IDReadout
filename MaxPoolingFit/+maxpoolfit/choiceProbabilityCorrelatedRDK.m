function [pRDK1,details] = choiceProbabilityCorrelatedRDK( ...
    cache,trialIdx,adaptation,noiseSD,rdkCorrelation,cfg,varargin)
%CHOICEPROBABILITYCORRELATEDRDK Integrate paired cross-RDK noise draws.
%
% Each RDK receives a directionally correlated Gaussian field with marginal
% SD noiseSD. rdkCorrelation controls the fraction of variance shared by the
% two RDKs. Corresponding draws are compared as pairs; all-cross-pairs
% integration would destroy the requested cross-RDK correlation.

p=inputParser;
p.addParameter('NDraws',cfg.nNoiseDraws, ...
  @(x)isnumeric(x)&&isscalar(x)&&x>=2&&mod(x,1)==0);
p.addParameter('ReturnWinners',false,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
rho=double(rdkCorrelation);
assert(isscalar(rho)&&isfinite(rho)&&rho>=0&&rho<=1, ...
  'RDK noise correlation must be a scalar between zero and one.');

nDraws=p.Results.NDraws;
zCommon=maxpoolfit.noiseDraws(cfg,nDraws,0);
zIndependent1=maxpoolfit.noiseDraws(cfg,nDraws,1000003);
zIndependent2=maxpoolfit.noiseDraws(cfg,nDraws,2000006);
z1=sqrt(rho).*zCommon+sqrt(1-rho).*zIndependent1;
z2=sqrt(rho).*zCommon+sqrt(1-rho).*zIndependent2;
[B,A,P]=maxpoolfit.deterministicResponse(cache,trialIdx,cfg);
pairedCfg=cfg; pairedCfg.useAllCrossPairs=false;
[pRDK1,details]=maxpoolfit.choiceProbabilityFromResponses( ...
  B,A,P,adaptation,noiseSD,pairedCfg,z1,z2,p.Results.ReturnWinners);
end
