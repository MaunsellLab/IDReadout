function sessionHeader = makeIDLSessionHeader(header,trials)
% makeIDLSessionHeader  Build lean session-level metadata for IDL analyses.
%
% The compact IDL session analyses inherit this structure as a top-level
% variable. Trial-dependent quantities, including the coherence steps that
% were actually used, are deliberately derived from the recorded trials.

sessionHeader = struct();
sessionHeader.version = 2;
sessionHeader.sourceDomain = 'IDL';
subjectNames = {'Meetz', 'Neesha', 'Test'};
isSubject = cellfun(@(s) contains(header.fileName, s), subjectNames);
assert(sum(isSubject) == 1, 'makeIDLSessionHeader:UnknownAnimal', 'Could not identify exactly one animal from %s.', ...
        header.fileName);
sessionHeader.animal = subjectNames{isSubject};
sessionHeader.date = localDataValue(header.date);
sessionHeader.numberOfTrials = localDataValue(header.numberOfTrials);
sessionHeader.frameRateHz = localDataValue(header.frameRateHz);
sessionHeader.preStepMS = localDataValue(header.preStepMS);
sessionHeader.stepMS = localDataValue(header.stepMS);



[~, sessionHeader.fileName] = fileparts(header.fileName);
sessionHeader.createdBy = mfilename;
sessionHeader.createdDate = datetime('now');
sessionHeader.metadataContract = 'lean session-level IDL analysis metadata';
if nargin>=2
  sessionHeader = idlSetAllocationHeader(sessionHeader,trials);
else
  assert(isfield(header,'blockStatus'), ...
    'makeIDLSessionHeader:MissingBlockStatus','Header has no allocation configuration.');
  info = idlAllocationConfiguration(header.blockStatus);
  names = fieldnames(info);
  for k=1:numel(names), sessionHeader.(names{k})=info.(names{k}); end
  sessionHeader.splitUnsplitRatios = info.splitUnsplitRatio;
  % Only a scan of all trials may mark a header as fully validated.
  sessionHeader.idlAllocationVersion = 0;
end
end

%% ------------------------------------------------------------------------
function value = localDataValue(value)

while isstruct(value) && isfield(value, 'data')
  value = value.data(1);
end
if isnumeric(value) || islogical(value) || isstring(value)
  assert(isscalar(value), 'makeIDLSessionHeader:ExpectedScalar', ...
    'Required numeric, logical, and string session-header values must be scalar.');
else
  assert(ischar(value), 'makeIDLSessionHeader:UnexpectedValueType', ...
    'Required session-header values must be numeric, logical, char, or string.');
end
end
