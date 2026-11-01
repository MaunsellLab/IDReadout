function H = idlRefreshSessionHeader(sourcePath)
% Backfill selection metadata in already converted FullSession MAT files.
S = load(sourcePath,'sessionHeader');
assert(isfield(S,'sessionHeader'), 'idlRefreshSessionHeader:MissingHeader', ...
  'FullSession lacks sessionHeader: %s',sourcePath);
H = S.sessionHeader;
if isfield(H,'idlAllocationVersion') && H.idlAllocationVersion==1 && ...
    all(isfield(H,{'twoAllocationMode','allocationStyle', ...
    'enabledAllocationPC','splitUnsplitRatio','splitUnsplitRatios'}))
  return
end
S = load(sourcePath,'trials');
assert(isfield(S,'trials'), 'idlRefreshSessionHeader:MissingTrials', ...
  'Cannot derive allocation metadata without recorded trials: %s',sourcePath);
sessionHeader = idlSetAllocationHeader(H,S.trials);
save(sourcePath,'sessionHeader','-append');
H = sessionHeader;
end
