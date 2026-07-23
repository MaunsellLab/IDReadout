function runID=formatDOGRunID(prefix,sigmaCenter,sigmaSurround,magnitude)
% formatDOGRunID  Stable filename-safe ID for one DOG grid member.
prefix=string(prefix);
if ~isscalar(prefix)||strlength(prefix)==0
  error('formatDOGRunID:BadPrefix','prefix must be a nonempty scalar string.');
end
if ~positiveScalar(sigmaCenter)||~nonnegativeScalar(magnitude)|| ...
    ~(isnan(sigmaSurround)||positiveScalar(sigmaSurround))
  error('formatDOGRunID:BadParameter','DOG parameters are invalid.');
end
if magnitude==0
  token=prefix+"_sigmaC"+numberToken(sigmaCenter)+"_baseline";
else
  if sigmaSurround<=sigmaCenter
    error('formatDOGRunID:BadSurround','Surround sigma must exceed center sigma.');
  end
  token=prefix+"_sigmaC"+numberToken(sigmaCenter)+ ...
    "_sigmaS"+numberToken(sigmaSurround)+"_a"+numberToken(magnitude);
end
if strlength(token)>80
  error('formatDOGRunID:TooLong','Generated RunID exceeds 80 characters.');
end
runID=token;
end

function token=numberToken(x)
token=replace(string(sprintf('%.9g',double(x))),["-","."],["m","p"]);
end
function tf=positiveScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;
end
function tf=nonnegativeScalar(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0;
end
