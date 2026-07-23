function summary=analyzeMTNegativeAsymptoteInvariance(varargin)
% analyzeMTNegativeAsymptoteInvariance  Test a constant negative readout tail.
%
% In the mean-subtracted MT model, adding a direction-independent term to
% every rotated readout is orthogonal to every stimulus and upstream-noise
% template. With unit center normalization,
%   W_b=(1+b)G-b,
% candidate signals and noise are therefore only scaled by 1+b. Positive
% scaling leaves rectification, winners, hard maxima, and choices unchanged.

p=inputParser;
addParameter(p,'NegativeAsymptoteMagnitude',[0 .02 .05 .1 .2 .4], ...
  @(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=0));
addParameter(p,'SigmaMTDeg',37.5,@positiveScalar);
addParameter(p,'SigmaReadoutDeg',5,@positiveScalar);
addParameter(p,'CandidateSpacingDeg',2,@positiveScalar);
addParameter(p,'InputSpacingDeg',1,@positiveScalar);
addParameter(p,'Seed',1729,@validSeed);
parse(p,varargin{:});opts=p.Results;

mt=makeMTReadoutForwardModel('sigmaMTDeg',opts.SigmaMTDeg);
candidateDirections=(0:opts.CandidateSpacingDeg:360-opts.CandidateSpacingDeg)';
inputDirections=(0:opts.InputSpacingDeg:360-opts.InputSpacingDeg)';
T=mtPopulationTemplate(inputDirections,mt);
G=makeGaussianReadoutBank(candidateDirections,opts.SigmaReadoutDeg,mt);
K0=T*double(G.weightsPhi)';C0=K0'*K0;

stream=RandStream('mt19937ar','Seed',opts.Seed);nTrial=200;
uSignalChange=randn(stream,nTrial,numel(inputDirections));
uSignalNoChange=randn(stream,nTrial,numel(inputDirections));
uNoiseChange=randn(stream,nTrial,numel(inputDirections));
uNoiseNoChange=randn(stream,nTrial,numel(inputDirections));
signalChange=uSignalChange*K0;signalNoChange=uSignalNoChange*K0;
noiseChange=uNoiseChange*K0;noiseNoChange=uNoiseNoChange*K0;
baseChoice=max(max(signalChange+noiseChange,0),[],2)> ...
  max(max(signalNoChange+noiseNoChange,0),[],2);

b=double(opts.NegativeAsymptoteMagnitude(:));n=numel(b);
expectedEvidenceScale=1+b;maxAbsTemplateSum=repmat(max(abs(sum(T,2))),n,1);
maxAbsSignalScaleResidual=nan(n,1);maxAbsCovarianceScaleResidual=nan(n,1);
choiceAgreement=nan(n,1);minimumWeight=nan(n,1);maximumWeight=nan(n,1);
for k=1:n
  B=makeGaussianReadoutBank(candidateDirections,opts.SigmaReadoutDeg,mt, ...
    'NegativeAsymptoteMagnitude',b(k));
  K=T*double(B.weightsPhi)';C=K'*K;scale=expectedEvidenceScale(k);
  maxAbsSignalScaleResidual(k)=max(abs(K-scale.*K0),[],'all');
  maxAbsCovarianceScaleResidual(k)=max(abs(C-scale.^2.*C0),[],'all');
  choice=max(max(uSignalChange*K+uNoiseChange*K,0),[],2)> ...
    max(max(uSignalNoChange*K+uNoiseNoChange*K,0),[],2);
  choiceAgreement(k)=mean(choice==baseChoice);
  minimumWeight(k)=min(B.weightsPhi,[],'all');
  maximumWeight(k)=max(B.weightsPhi,[],'all');
end

summary=table(b,minimumWeight,maximumWeight,expectedEvidenceScale, ...
  maxAbsTemplateSum,maxAbsSignalScaleResidual, ...
  maxAbsCovarianceScaleResidual,choiceAgreement, ...
  'VariableNames',{'negativeAsymptoteMagnitude','minimumWeight', ...
  'maximumWeight','expectedEvidenceScale','maxAbsTemplateSum', ...
  'maxAbsSignalScaleResidual','maxAbsCovarianceScaleResidual', ...
  'choiceAgreement'});
disp('MT constant-negative-asymptote invariance test');disp(summary);
end

function tf=positiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=validSeed(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=2^32-1&&fix(x)==x;
end
