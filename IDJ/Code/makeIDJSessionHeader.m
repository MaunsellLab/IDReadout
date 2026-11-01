function sessionHeader = makeIDJSessionHeader(header)
% makeIDJSessionHeader  Build lean session-level metadata for IDJ analyses.
%
% The compact IDJ session analyses inherit this structure as a top-level
% variable. Trial-dependent quantities, including the coherence steps that
% were actually used, are deliberately derived from the recorded trials.

sessionHeader = struct();
sessionHeader.version = 1;
sessionHeader.sourceDomain = 'IDJ';
subjectNames = {'Meetz', 'Neesha', 'Test'};
isSubject = cellfun(@(s) contains(header.fileName, s), subjectNames);
assert(sum(isSubject) == 1, 'makeIDJSessionHeader:UnknownAnimal', 'Could not identify exactly one animal from %s.', ...
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
sessionHeader.metadataContract = 'lean session-level IDJ analysis metadata';
end

%% ------------------------------------------------------------------------
function value = localDataValue(value)

while isstruct(value) && isfield(value, 'data')
  value = value.data(1);
end
if isnumeric(value) || islogical(value) || isstring(value)
  assert(isscalar(value), 'makeIDJSessionHeader:ExpectedScalar', ...
    'Required numeric, logical, and string session-header values must be scalar.');
else
  assert(ischar(value), 'makeIDJSessionHeader:UnexpectedValueType', ...
    'Required session-header values must be numeric, logical, char, or string.');
end
end
