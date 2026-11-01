function cache = loadCache(path)
%LOADCACHE Load and validate a MaxPoolingTrialCache-v2 MAT file.
S = load(path, 'cache');
assert(isfield(S,'cache'), 'File does not contain a variable named cache.');
cache = S.cache;
maxpoolfit.validateCache(cache);
end
