function tag = idlAllocationTag(twoAllocationMode)
% Preserve original output names; isolate the new analysis outputs.
tag = '';
if twoAllocationMode, tag = '_TwoAlloc'; end
end
