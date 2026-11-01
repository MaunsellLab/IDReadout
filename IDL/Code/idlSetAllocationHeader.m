function H = idlSetAllocationHeader(H,trials)
% Derive header selection metadata from every recorded trial, including aborts.
assert(iscell(trials) && ~isempty(trials), ...
  'idlSetAllocationHeader:Trials','Expected a nonempty trial cell array.');
ratios = nan(numel(trials),1);
for k=1:numel(trials)
  tr = trials{k};
  assert(isfield(tr,'blockStatus'), 'idlSetAllocationHeader:MissingBlockStatus', ...
    'Trial %d has no blockStatus configuration.',k);
  info = idlAllocationConfiguration(tr.blockStatus);
  if k==1
    first = info;
  else
    assert(info.twoAllocationMode==first.twoAllocationMode, ...
      'idlSetAllocationHeader:MixedStyles', ...
      'A file cannot mix original and two-allocation analysis styles (trial %d).',k);
  end
  ratios(k) = info.splitUnsplitRatio;
  if isfield(tr,'trial')
    D = tr.trial.data; B = tr.blockStatus.data(1);
    index = double(D.allocIndex)+1;
    assert(isscalar(index) && index==round(index) && index>=1 && ...
      index<=double(B.numAllocs) && logical(B.allocEnabled(index)) && ...
      abs(double(D.alloc0)-double(B.alloc0(index)))<1e-6 && ...
      abs(double(D.alloc1)-double(B.alloc1(index)))<1e-6, ...
      'idlSetAllocationHeader:TrialConfiguration', ...
      'Trial %d disagrees with its enabled allocation slot.',k);
  end
end
if ~isfield(H,'version'), H.version=1; end
H.version = max(2,double(H.version));
H.sourceDomain = 'IDL';
H.idlAllocationVersion = 1;
H.twoAllocationMode = first.twoAllocationMode;
H.allocationStyle = first.allocationStyle;
H.enabledAllocationPC = first.enabledAllocationPC;
H.splitUnsplitRatios = unique(ratios)';
H.splitUnsplitRatio = NaN;
if numel(H.splitUnsplitRatios)==1, H.splitUnsplitRatio=H.splitUnsplitRatios; end
end
