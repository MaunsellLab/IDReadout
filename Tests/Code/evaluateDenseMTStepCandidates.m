function result=evaluateDenseMTStepCandidates(componentDirectionsDeg, ...
  componentCoherence,bank,mtModel,varargin)
% evaluateDenseMTStepCandidates  Evidence for a fixed dense readout bank.

p=inputParser;
addParameter(p,'CandidateNoise',[],@(x) isempty(x)||isnumeric(x));
addParameter(p,'RectifyCandidates',true,@(x) islogical(x)&&isscalar(x));
addParameter(p,'MTTemplateMode',"signed",@validTemplateMode);
parse(p,varargin{:}); opts=p.Results;

coherence=double(componentCoherence);
if isvector(coherence),coherence=coherence(:)';end
templates=mtPopulationTemplate(componentDirectionsDeg,mtModel, ...
  'MTTemplateMode',opts.MTTemplateMode);
componentToCandidate=templates*double(bank.weightsPhi)';
evidence=coherence*componentToCandidate;
if ~isempty(opts.CandidateNoise)
  if ~isequal(size(opts.CandidateNoise),size(evidence))
    error('evaluateDenseMTStepCandidates:NoiseSizeMismatch', ...
      'CandidateNoise must match observations by candidates.');
  end
  evidence=evidence+double(opts.CandidateNoise);
end
activation=evidence;
if opts.RectifyCandidates,activation=max(activation,0);end
[maximum,winner]=max(activation,[],2);

result=struct('evidence',evidence,'activation',activation, ...
  'maximum',maximum,'winner',winner, ...
  'candidateDirectionsDeg',bank.candidateDirectionsDeg, ...
  'componentToCandidate',componentToCandidate, ...
  'mtTemplateMode',lower(string(opts.MTTemplateMode)), ...
  'rectifyCandidates',opts.RectifyCandidates);
end

function tf=validTemplateMode(x)
tf=(ischar(x)||isstring(x))&&isscalar(string(x))&& ...
  any(strcmpi(string(x),["signed","rectified"]));
end
