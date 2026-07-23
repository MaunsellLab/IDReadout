function row=matchSyntheticTrialRows(sourceSession,sourceTrial,sourceOffset, ...
  syntheticSession,syntheticTrial,syntheticOffset)
% matchSyntheticTrialRows  Exact one-to-one trial-key matching.

% Offset vectors may be empty when session and trial identify rows uniquely.

sourceSession=double(sourceSession(:));sourceTrial=double(sourceTrial(:));
syntheticSession=double(syntheticSession(:));syntheticTrial=double(syntheticTrial(:));
if isempty(sourceOffset)
  sourceOffset=zeros(size(sourceSession));syntheticOffset=zeros(size(syntheticSession));
else
  sourceOffset=double(sourceOffset(:));syntheticOffset=double(syntheticOffset(:));
end
if numel(sourceSession)~=numel(sourceTrial)||numel(sourceSession)~=numel(sourceOffset)
  error('matchSyntheticTrialRows:SourceSizeMismatch', ...
    'Source key columns must have equal lengths.');
end
if numel(syntheticSession)~=numel(syntheticTrial)|| ...
    numel(syntheticSession)~=numel(syntheticOffset)
  error('matchSyntheticTrialRows:SyntheticSizeMismatch', ...
    'Synthetic key columns must have equal lengths.');
end
if any(~isfinite([sourceSession;sourceTrial;sourceOffset]))|| ...
    any(~isfinite([syntheticSession;syntheticTrial;syntheticOffset]))
  error('matchSyntheticTrialRows:NonfiniteKey','Trial keys must be finite.');
end
sourceKey=makeKey(sourceSession,sourceTrial,sourceOffset);
syntheticKey=makeKey(syntheticSession,syntheticTrial,syntheticOffset);
if numel(unique(sourceKey))~=numel(sourceKey)
  error('matchSyntheticTrialRows:DuplicateSourceKey', ...
    'Source trial keys are not unique.');
end
if numel(unique(syntheticKey))~=numel(syntheticKey)
  error('matchSyntheticTrialRows:DuplicateSyntheticKey', ...
    'Synthetic trial keys are not unique.');
end
[found,row]=ismember(syntheticKey,sourceKey);
if ~all(found)
  error('matchSyntheticTrialRows:MissingTrial', ...
    '%d synthetic trials were absent from the source table.',sum(~found));
end
row=double(row(:));
end

function key=makeKey(session,trial,offset)
key=compose('%.17g|%.17g|%.17g',session,trial,offset);
end
