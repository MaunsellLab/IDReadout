function [nBytes,fileTable] = simulationRunSize(context)
% simulationRunSize  Return the total bytes and inventory for one run.

context = validateAnalysisContext(context);
if context.Mode ~= "synthetic"
  error("IDReadout:SyntheticContextRequired", ...
    "simulationRunSize requires a synthetic context.");
end
runRoot = string(domainFolder(mfilename('fullpath'),'Tests',context));
if ~isfolder(runRoot)
  nBytes = 0;
  fileTable = table(strings(0,1),zeros(0,1), ...
    'VariableNames',{'relativePath','bytes'});
  return
end
files = dir(fullfile(runRoot,'**','*'));
files = files(~[files.isdir]);
relativePath = strings(numel(files),1);
bytes = zeros(numel(files),1);
prefix = runRoot + filesep;
for k = 1:numel(files)
  fullPath = string(fullfile(files(k).folder,files(k).name));
  relativePath(k) = extractAfter(fullPath,prefix);
  bytes(k) = files(k).bytes;
end
fileTable = table(relativePath,bytes);
nBytes = sum(bytes);
end
