function [exclude, reasons] = applyExperimentalValidityChecks(metadata)
% applyExperimentalValidityChecks  Apply current session-level validity rules.
%
% This wrapper gives the call site a name that describes the purpose of the
% check while preserving excludeFile.m as the current implementation.

[exclude, reason] = excludeFile(metadata);
if exclude
  reasons = {reason};
else
  reasons = {};
end
end
