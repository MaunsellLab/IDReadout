function bank = makeGaussianReadoutBank(candidateDirectionsDeg,sigmaReadoutDeg,mtModel,varargin)
% makeGaussianReadoutBank  Rotated Gaussian MT readouts.

% Each row of bank.weightsPhi is centered on one candidate direction and
% has unit peak. NegativeAsymptoteMagnitude b defines
%   W=(1+b)*Gaussian-b,
% so W(0)=1 and the far-direction asymptote is -b. The default b=0 is the
% original nonnegative Gaussian readout.

p=inputParser;
addParameter(p,'NegativeAsymptoteMagnitude',0,@isNonnegativeScalar);
addParameter(p,'SurroundMagnitude',0,@isNonnegativeScalar);
addParameter(p,'SigmaSurroundDeg',nan,@isOptionalPositiveScalar);
parse(p,varargin{:});opts=p.Results;

if ~isnumeric(candidateDirectionsDeg) || ~isvector(candidateDirectionsDeg) || ...
    any(~isfinite(candidateDirectionsDeg))
  error('makeGaussianReadoutBank:BadDirections', ...
    'candidateDirectionsDeg must be a finite numeric vector.');
end
if ~isnumeric(sigmaReadoutDeg) || ~isscalar(sigmaReadoutDeg) || ...
    ~isfinite(sigmaReadoutDeg) || sigmaReadoutDeg<=0
  error('makeGaussianReadoutBank:BadSigma', ...
    'sigmaReadoutDeg must be a finite positive scalar.');
end
if ~isstruct(mtModel) || ~isfield(mtModel,'phiDeg')
  error('makeGaussianReadoutBank:BadMTModel','mtModel must contain phiDeg.');
end

directions = mod(double(candidateDirectionsDeg(:))+180,360)-180;
phiDeg = double(mtModel.phiDeg(:)');
distanceDeg = mod(phiDeg-directions+180,360)-180;
gaussianPhi = exp(-(distanceDeg.^2)./(2*double(sigmaReadoutDeg).^2));
b=double(opts.NegativeAsymptoteMagnitude);
s=double(opts.SurroundMagnitude);
sigmaSurround=double(opts.SigmaSurroundDeg);
if s>0
  if ~isfinite(sigmaSurround)||sigmaSurround<=double(sigmaReadoutDeg)
    error('makeGaussianReadoutBank:BadSurroundSigma', ...
      'SigmaSurroundDeg must exceed sigmaReadoutDeg when the surround is nonzero.');
  end
  surroundPhi=exp(-(distanceDeg.^2)./(2*sigmaSurround.^2));
else
  surroundPhi=zeros(size(gaussianPhi));
end
weightsPhi=(1+b+s).*gaussianPhi-s.*surroundPhi-b;

bank = struct();
bank.candidateDirectionsDeg = directions;
bank.sigmaReadoutDeg = double(sigmaReadoutDeg);
bank.negativeAsymptoteMagnitude=b;
bank.surroundMagnitude=s;
bank.sigmaSurroundDeg=sigmaSurround;
bank.phiDeg = phiDeg;
bank.weightsPhi = weightsPhi;
bank.definition = sprintf([ ...
  'Unit-center Gaussian readout: (1+b+s)Gcenter-sGsurround-b; ' ...
  'negative asymptote magnitude %.9g; surround magnitude %.9g; ' ...
  'surround sigma %.9g deg'],b,s,sigmaSurround);
end

function tf=isNonnegativeScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0;
end
function tf=isOptionalPositiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&(isnan(x)||(isfinite(x)&&x>0));
end
