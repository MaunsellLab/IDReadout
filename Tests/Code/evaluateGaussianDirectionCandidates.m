function result = evaluateGaussianDirectionCandidates( ...
  populationResponse,bank,varargin)
% evaluateGaussianDirectionCandidates  Read and optionally pool MT responses.
%
% Name-value options:
%   RectifyCandidates  max(evidence,0) before pooling (default false)
%   PoolingMode        'none', 'max', 'sum', or 'pnorm' (default 'none')
%   P                  p for p-norm pooling (default 2)

p=inputParser;
addParameter(p,'RectifyCandidates',false,@(x) islogical(x)&&isscalar(x));
addParameter(p,'PoolingMode','none',@(x) ischar(x)||isstring(x));
addParameter(p,'P',2,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1);
parse(p,varargin{:});
opts=p.Results;

populationResponse=double(populationResponse);
if isvector(populationResponse), populationResponse=populationResponse(:)'; end
if size(populationResponse,2)~=size(bank.weightsPhi,2)
  error('evaluateGaussianDirectionCandidates:GridMismatch', ...
    'Population response and readout bank use different phi grids.');
end

evidence=populationResponse*double(bank.weightsPhi)';
activation=evidence;
if opts.RectifyCandidates, activation=max(activation,0); end

mode=lower(string(opts.PoolingMode));
switch mode
  case "none"
    pooled=[]; winner=[];
  case "max"
    [pooled,winner]=max(activation,[],2);
  case "sum"
    pooled=sum(activation,2); winner=[];
  case "pnorm"
    if any(activation<0,'all')
      error('evaluateGaussianDirectionCandidates:SignedPNorm', ...
        'P-norm pooling requires nonnegative candidate activations.');
    end
    pooled=sum(activation.^opts.P,2).^(1/opts.P); winner=[];
  otherwise
    error('evaluateGaussianDirectionCandidates:BadPoolingMode', ...
      'Unknown pooling mode: %s',mode);
end

result=struct();
result.evidence=evidence;
result.activation=activation;
result.pooled=pooled;
result.winner=winner;
result.rectifyCandidates=opts.RectifyCandidates;
result.poolingMode=mode;
result.p=opts.P;
end
