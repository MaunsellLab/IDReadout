function [summaries,notes] = findLocalSummaries(sourceFolders,animal)
%FINDLOCALSUMMARIES Find saved, animal-specific, native IDS gain summaries.
summaries={}; notes=strings(0,1); paths=strings(0,1);
for k=1:numel(sourceFolders)
  folder=char(sourceFolders{k}); [parent,leaf]=fileparts(folder);
  if strcmpi(leaf,'ProbeSessions'), base=parent; else, base=folder; end
  D=dir(fullfile(base,'AcrossSessionSummaries','*.mat'));
  for q=1:numel(D), paths(end+1,1)=string(fullfile(D(q).folder,D(q).name)); end %#ok<AGROW>
end
paths=unique(paths,'stable');
for k=1:numel(paths)
  names=who('-file',paths(k));
  if ~ismember('acrossSummary',names), continue; end
  A=load(paths(k),'acrossSummary'); S=A.acrossSummary;
  if ~isfield(S,'animal') || ~strcmpi(string(S.animal),animal)
    notes(end+1,1)="Skipped summary with different/unspecified animal: "+paths(k); %#ok<AGROW>
    continue
  end
  if ~isfield(S,'sessionRecords') || ...
      ~ismember('sourceDomain',S.sessionRecords.Properties.VariableNames)
    notes(end+1,1)="Skipped summary without explicit sourceDomain records: "+paths(k); %#ok<AGROW>
    continue
  end
  if ~all(strcmpi(string(S.sessionRecords.sourceDomain),'IDS'))
    notes(end+1,1)="Skipped summary containing non-IDS sessions: "+paths(k); %#ok<AGROW>
    continue
  end
  summaries{end+1,1}=S; %#ok<AGROW>
end
if isempty(summaries)
  notes(end+1,1)="Run plotIDSSummary('Animal','Meetz','IncludeIDQ',false), then recover/export gains.";
end
end
