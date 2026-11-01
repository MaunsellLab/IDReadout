function manifest = discoverProbeSessions(projectRoot, varargin)
%DISCOVERPROBESESSIONS Select and inventory experimental ProbeSessions.
%
% Candidate files are passed through selectAnalysisFiles with
% LoadMetadata=false. Thus its centralized hard exclusions are applied,
% while unchanged MAT files remain unopened during routine cache scans.

p=inputParser;
p.addParameter('Domains',["IDR" "IDQ" "IDS"],@(x)isstring(x)||iscellstr(x));
p.addParameter('Animals',["Meetz" "Neesha"],@(x)isstring(x)||iscellstr(x));
p.addParameter('FileSelectionArgs',{},@iscell);
p.parse(varargin{:});
domains=upper(string(p.Results.Domains)); animals=string(p.Results.Animals);
selectionArgs=p.Results.FileSelectionArgs;
validateSelectionArgs(selectionArgs);
projectRoot=string(projectRoot);
assert(isfolder(projectRoot),'Project root not found: %s',projectRoot);
pathCleanup=ensureSelectorAvailable(projectRoot); %#ok<NASGU>

paths=strings(0,1); domainCol=strings(0,1); animalCol=strings(0,1);
bytes=zeros(0,1); modifiedDatenum=zeros(0,1);
for domain=domains(:)'
  [folders,filePattern]=probeFolders(projectRoot,domain);
  if isempty(folders), continue; end
  selected=selectAnalysisFiles(cellstr(folders), ...
    'FilePattern',filePattern,'LoadMetadata',false, ...
    'Bin180Into179',false,'IncludeIDQ',false,selectionArgs{:});
  selected=string(selected(:));
  for k=1:numel(selected)
    filePath=selected(k);
    hit=arrayfun(@(animal)contains(filePath,animal),animals);
    if ~any(hit), continue; end
    assert(sum(hit)==1,'Cannot uniquely identify animal for %s.',filePath);
    D=dir(char(filePath));
    assert(isscalar(D)&&~D.isdir,'Selected source is missing or ambiguous: %s',filePath);
    paths(end+1,1)=filePath; %#ok<AGROW>
    domainCol(end+1,1)=domain; %#ok<AGROW>
    animalCol(end+1,1)=animals(hit); %#ok<AGROW>
    bytes(end+1,1)=double(D.bytes); %#ok<AGROW>
    modifiedDatenum(end+1,1)=double(D.datenum); %#ok<AGROW>
  end
end
manifest=table(paths,domainCol,animalCol,bytes,modifiedDatenum, ...
  'VariableNames',{'sourceFile','domain','animal','bytes','modifiedDatenum'});
manifest=sortrows(manifest,{'domain','animal','sourceFile'});
assert(numel(unique(manifest.sourceFile))==height(manifest), ...
  'Directory scan returned duplicate source paths.');
end

function [folders,filePattern]=probeFolders(projectRoot,domain)
switch domain
  case {"IDR","IDS"}
    D=dir(char(fullfile(projectRoot,domain,'Data','Probe*')));
    D=D([D.isdir]); folders=strings(0,1);
    for k=1:numel(D)
      candidate=string(fullfile(D(k).folder,D(k).name,'ProbeSessions'));
      if isfolder(candidate), folders(end+1,1)=candidate; end %#ok<AGROW>
    end
    filePattern='*.mat';
  case "IDQ"
    D=dir(char(fullfile(projectRoot,'IDQ','Data','Probe*')));
    D=D([D.isdir]); folders=strings(0,1);
    for k=1:numel(D)
      folders(end+1,1)=string(fullfile(D(k).folder,D(k).name)); %#ok<AGROW>
    end
    filePattern='ProbeSession*.mat';
  otherwise
    error('Unsupported domain %s.',domain);
end
end

function cleanup=ensureSelectorAvailable(projectRoot)
candidate=char(fullfile(projectRoot,'selectAnalysisFiles.m'));
if isfile(candidate)
  resolved=which('selectAnalysisFiles');
  if strcmp(resolved,candidate), cleanup=[]; return; end
  oldPath=path; addpath(char(projectRoot),'-begin');
  cleanup=onCleanup(@()path(oldPath));
  assert(strcmp(which('selectAnalysisFiles'),candidate), ...
    'Could not activate the project-root selectAnalysisFiles.m: %s',candidate);
  return
end
assert(exist('selectAnalysisFiles','file')==2, ...
  ['selectAnalysisFiles.m was not found on the MATLAB path or in the ' ...
   'project root: %s'],projectRoot);
cleanup=[];
end

function validateSelectionArgs(args)
assert(mod(numel(args),2)==0,'FileSelectionArgs must contain name/value pairs.');
allowed=["HardExcludeFileNames" "HardExcludeNameStrings" "ReportExcluded"];
for k=1:2:numel(args)
  assert(ischar(args{k})||isstring(args{k}),'FileSelectionArgs names must be text.');
  assert(any(strcmpi(string(args{k}),allowed)), ...
    'Unsupported filename-only FileSelectionArgs option: %s',string(args{k}));
end
end
