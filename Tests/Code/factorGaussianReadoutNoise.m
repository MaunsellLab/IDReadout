function noiseModel = factorGaussianReadoutNoise(bank,varargin)
% factorGaussianReadoutNoise  Correlated readout noise from white MT noise.
%
% If eta_phi has independent unit-variance elements, readout noise is
% eta_phi*W'. Its covariance is W*W'. The retained eigensystem provides an
% efficient exact-to-tolerance sampler without storing MT-neuron noise.

p=inputParser;
addParameter(p,'RelativeEigenvalueTolerance',1e-10, ...
  @(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1);
parse(p,varargin{:}); opts=p.Results;

W=double(bank.weightsPhi);
C=(W*W'); C=(C+C')/2;
[V,D]=eig(C,'vector');
[D,order]=sort(real(D),'descend'); V=real(V(:,order));
tolerance=max(D)*opts.RelativeEigenvalueTolerance;
keep=D>max(tolerance,0);
factor=V(:,keep).*sqrt(D(keep))';

noiseModel=struct();
noiseModel.covariance=C;
noiseModel.factor=factor;
noiseModel.rank=sum(keep);
noiseModel.relativeEigenvalueTolerance=opts.RelativeEigenvalueTolerance;
noiseModel.candidateDirectionsDeg=bank.candidateDirectionsDeg;
noiseModel.definition=[ ...
  "White homoscedastic Gaussian variability over the MT preferred-" + ...
  "direction population, projected through the Gaussian readout bank"];
end
