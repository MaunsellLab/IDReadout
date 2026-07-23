function runID=formatReadoutWidthRunID(prefix,sigmaReadoutDeg)
% formatReadoutWidthRunID  Stable filesystem-safe width-run identifier.
mustBeTextScalar(prefix);prefix=string(prefix);
if strlength(prefix)==0
  error('formatReadoutWidthRunID:EmptyPrefix','RunIDPrefix cannot be empty.');
end
if ~(isnumeric(sigmaReadoutDeg)&&isscalar(sigmaReadoutDeg)&& ...
    isfinite(sigmaReadoutDeg)&&sigmaReadoutDeg>0)
  error('formatReadoutWidthRunID:BadWidth', ...
    'sigmaReadoutDeg must be a finite positive scalar.');
end
token=string(sprintf('%.12g',double(sigmaReadoutDeg)));
token=replace(token,'.','p');token=replace(token,'-','m');token=replace(token,'+','');
runID=prefix+"_sigmaR"+token;
if strlength(runID)>80
  error('formatReadoutWidthRunID:RunIDTooLong', ...
    'Width-specific RunID exceeds the 80-character context limit.');
end
end
