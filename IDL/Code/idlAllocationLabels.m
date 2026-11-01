function labels = idlAllocationLabels(allocationPC)
% Format one allocation per label, independently of input row/column shape.
labels = arrayfun(@(a) sprintf('%g:%g',a,100-a), ...
  allocationPC,'UniformOutput',false);
end
