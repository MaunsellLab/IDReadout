function noiseModel=factorUpstreamDirectionalNoise( ...
  inputDirectionsDeg,bank,mtModel,varargin)
% factorUpstreamDirectionalNoise  Pass directional input noise through MT.
%
% A unit upstream perturbation u(direction) produces MT response u*T and
% candidate evidence u*T*W'. No independent MT-response or readout noise is
% added. The eigensystem is an efficient sampler of the resulting candidate
% covariance.

p=inputParser;
addParameter(p,'RelativeEigenvalueTolerance',1e-10, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
parse(p,varargin{:});opts=p.Results;

T=mtPopulationTemplate(inputDirectionsDeg,mtModel); % input direction x MT
inputToCandidate=T*double(bank.weightsPhi)';
candidateCovariance=inputToCandidate'*inputToCandidate;
candidateCovariance=(candidateCovariance+candidateCovariance')/2;
[V,D]=eig(candidateCovariance,'vector');
[D,order]=sort(real(D),'descend');V=real(V(:,order));
tolerance=max(D)*opts.RelativeEigenvalueTolerance;
keep=D>max(tolerance,0);
factor=V(:,keep).*sqrt(D(keep))';

mtCovariance=T'*T;
mtSD=sqrt(diag(mtCovariance));
mtCorrelation=mtCovariance./(mtSD*mtSD');

noiseModel=struct();
noiseModel.inputDirectionsDeg=double(inputDirectionsDeg(:));
noiseModel.inputToMT=T;
noiseModel.inputToCandidate=inputToCandidate;
noiseModel.mtCovariance=mtCovariance;
noiseModel.mtCorrelation=mtCorrelation;
noiseModel.candidateCovariance=candidateCovariance;
noiseModel.factor=factor;
noiseModel.rank=sum(keep);
noiseModel.relativeEigenvalueTolerance=opts.RelativeEigenvalueTolerance;
noiseModel.definition=[ ...
  "Independent homoscedastic Gaussian perturbations over latent motion " + ...
  "input directions, passed through MT tuning and the supplied readouts; no " + ...
  "separate MT-response or readout noise"];
end
