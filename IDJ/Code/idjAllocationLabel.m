function label = idjAllocationLabel(majorPC)
% idjAllocationLabel  Label one nominal three-stream IDJ allocation.
if abs(majorPC-100/3)<1e-6
  label = '33.3:33.3:33.3';
else
  label = sprintf('%g:%g:%g',majorPC,(100-majorPC)/2,(100-majorPC)/2);
end
end
