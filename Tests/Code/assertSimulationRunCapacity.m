function [currentBytes,remainingBytes] = assertSimulationRunCapacity(context,additionalBytes)
% assertSimulationRunCapacity  Enforce the context's per-run disk quota.

context = validateAnalysisContext(context);
if context.Mode ~= "synthetic"
  error("IDReadout:SyntheticContextRequired", ...
    "Run-capacity checks require a synthetic context.");
end
if ~isnumeric(additionalBytes) || ~isscalar(additionalBytes) || ...
    ~isfinite(additionalBytes) || additionalBytes<0
  error("IDReadout:InvalidAdditionalBytes", ...
    "additionalBytes must be a finite nonnegative scalar.");
end
currentBytes = simulationRunSize(context);
remainingBytes = context.MaxRunBytes-currentBytes;
if currentBytes+additionalBytes > context.MaxRunBytes
  error("IDReadout:SyntheticRunQuotaExceeded", ...
    "Run %s would exceed its %.3f GB quota. Current size %.3f GB; " + ...
    "requested additional size %.3f GB.",context.RunID, ...
    context.MaxRunBytes/1e9,currentBytes/1e9,additionalBytes/1e9);
end
end
