function [files,info,folders] = selectLocalFiles(root,animal,reportExcluded)
%SELECTLOCALFILES Delegate metadata and exclusions to the project selector.
assert(exist('selectAnalysisFiles','file')==2, ...
  'Add the current project selectAnalysisFiles.m to the MATLAB path.');
root=char(root); assert(isfolder(root),'Data folder not found: %s',root);
if isfolder(fullfile(root,'Data')), base=fullfile(root,'Data'); else, base=root; end
[~,leaf]=fileparts(root);
if strcmpi(leaf,'ProbeSessions')
  folders={root};
elseif isfolder(fullfile(root,'ProbeSessions'))
  folders={fullfile(root,'ProbeSessions')};
else
  D=dir(fullfile(base,'Probe*')); D=D([D.isdir]); folders={};
  for k=1:numel(D)
    folder=fullfile(D(k).folder,D(k).name,'ProbeSessions');
    if isfolder(folder), folders{end+1,1}=folder; end %#ok<AGROW>
  end
  % A directly selected folder of converted raw sessions also works, when
  % it has the SourceDomain metadata required by the project selector.
  if isempty(folders), folders={base}; end
end
[files,info]=selectAnalysisFiles(folders, ...
  'Animal',char(animal),'SourceDomain','IDS','IncludeIDQ',false, ...
  'Bin180Into179',false,'FilePattern','*.mat','LoadMetadata',true, ...
  'ReportExcluded',reportExcluded);
files=cellstr(string(files(:)));
assert(~isempty(files),'No %s IDS files passed selectAnalysisFiles.',animal);
[files,order]=sort(files); info=info(order,:);
assert(numel(unique(string(files)))==numel(files),'Selector returned duplicate paths.');
end
