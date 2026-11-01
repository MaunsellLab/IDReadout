function s = pool(r,p)
%POOL Nonnegative rate p-norm, pooling down the first array dimension.
assert(isscalar(p) && p>=1 && ~isnan(p), 'p must be >=1 or Inf.');
assert(all(isfinite(r),'all') && all(r>=0,'all'), ...
  'Pool inputs must be finite, nonnegative firing responses.');
m=max(r,[],1);
if isinf(p), s=m; return; end
den=m; den(den==0)=1;
s=m.*sum((r./den).^p,1).^(1/p);
end
