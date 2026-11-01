function [items,notes] = discover(root,options)
%DISCOVER Select real source files; use their compact analyses when present.
% The project selector owns filename/path exclusions. Metadata are read here
% so legacy compact products with nested headers can also be exported.
assert(exist('selectAnalysisFiles','file')==2,'Add your current selectAnalysisFiles.m to the path.');
items=struct('source',{},'compact',{},'locationDomain',{},'domain',{}, ...
  'animal',{},'session',{},'offset',{},'randomBase',{});
notes=cell(0,3); roots={}; domains={};
requested=cellstr(string(options.Domains)); [~,rootLeaf]=fileparts(root);
for k=1:numel(requested)
  d=upper(requested{k}); candidate=fullfile(root,d);
  if isfolder(candidate)
    roots{end+1}=candidate; domains{end+1}=d; %#ok<AGROW>
  elseif strcmpi(d,'IDR')&&isfolder(fullfile(root,'Data'))&& ...
      ~any(strcmpi(rootLeaf,{'IDS','IDQ','IDJ'}))
    roots{end+1}=root; domains{end+1}=d; %#ok<AGROW>
  else
    [~,leaf]=fileparts(root);
    if strcmpi(leaf,d), roots{end+1}=root; domains{end+1}=d; end %#ok<AGROW>
  end
end
for r=1:numel(roots)
  dr=roots{r}; declared=domains{r}; base=fullfile(dr,'Data');
  if ~isfolder(base), base=dr; end
  folders={};
  if ~strcmpi(declared,'IDJ')
    D=dir(fullfile(base,'Probe*')); D=D([D.isdir]);
    for q=1:numel(D)
      p=fullfile(D(q).folder,D(q).name,'ProbeSessions');
      if isfolder(p), folders{end+1}=p;
      elseif strcmpi(declared,'IDQ')
        % Historical IDQ ProbeSession files live directly in Probe120.
        folders{end+1}=fullfile(D(q).folder,D(q).name);
      end %#ok<AGROW>
    end
  end
  if isempty(folders) && isfolder(fullfile(base,'FullSessions'))
    folders={fullfile(base,'FullSessions')};
  end
  if isempty(folders)
    notes(end+1,:)={declared,dr,'No ProbeSessions or FullSessions directory found'}; %#ok<AGROW>
    continue
  end
  for f=1:numel(folders)
    % Animal/domain filtering below uses explicit headers, not filenames.
    paths=selectAnalysisFiles(folders{f},'LoadMetadata',false, ...
      'IncludeIDQ',false,'Bin180Into179',false,'FilePattern','*.mat');
    for q=1:numel(paths)
      path=char(paths{q});
      try
        vars=who('-file',path);
        if ~any(ismember({'probeTrials','trials'},vars)), continue; end
        H=load(path,'sessionHeader'); assert(isfield(H,'sessionHeader'),'Missing sessionHeader');
        h=H.sessionHeader; p=struct();
        if ismember('sessionProbeHeader',vars)
          P=load(path,'sessionProbeHeader'); p=P.sessionProbeHeader;
        end
        animal=char(mtexport.get(h,{'animal','subject'},''));
        if isempty(animal), error('Source lacks animal metadata'); end
        if ~strcmpi(options.Animal,'All')&&~strcmpi(options.Animal,animal), continue; end
        domain=upper(char(mtexport.get(p,'sourceDomain',mtexport.get(h,'sourceDomain',declared))));
        if any(strcmpi(domain,{'IDReadout','IDReadout2'})), domain='IDR'; end
        if ~any(strcmpi(domain,requested)), continue; end
        randomBase=double(mtexport.get(p,'randomBaseDir',mtexport.get(h,'randomBaseDir',NaN)));
        if strcmpi(domain,'IDQ')&&isfinite(randomBase)&&~logical(randomBase)
          notes(end+1,:)={domain,path,'Excluded explicitly fixed-base-direction IDQ'}; %#ok<AGROW>
          continue
        end
        name=mtexport.get(p,'sourceSessionID',mtexport.get(p,'parentFileName',mtexport.get(h,'fileName',path)));
        sid=mtexport.sessionID(name);
        offset=double(mtexport.get(p,{'probeOffsetDeg','probeDirDeg'},NaN));
        [parent,leaf]=fileparts(folders{f}); [~,stem]=fileparts(path);
        compact='';
        if strcmpi(leaf,'ProbeSessions')
          candidate=fullfile(parent,'SessionAnalysis',[stem '_sessionAnalysis.mat']);
        elseif strcmpi(leaf,'FullSessions')
          candidate=fullfile(parent,'SessionAnalysis',[stem '_sessionAnalysis.mat']);
        else
          candidate=fullfile(folders{f},'SessionAnalysis',[stem '_sessionAnalysis.mat']);
        end
        if isfile(candidate), compact=candidate; end
        items(end+1)=struct('source',path,'compact',compact,'locationDomain',declared, ...
          'domain',domain,'animal',animal,'session',sid,'offset',offset,'randomBase',randomBase); %#ok<AGROW>
      catch ME
        notes(end+1,:)={declared,path,['Discovery error: ' ME.message]}; %#ok<AGROW>
      end
    end
  end
end
% Prefer role-normalized IDQ sources in the IDS tree to legacy equivalents.
if isempty(items), return; end
priority=ones(numel(items),1);
for k=1:numel(items)
  if strcmp(items(k).domain,'IDQ')&&strcmp(items(k).locationDomain,'IDS'), priority(k)=0; end
end
[~,order]=sort(priority,'ascend'); items=items(order);
seen=containers.Map('KeyType','char','ValueType','logical'); keep=true(size(items));
for k=1:numel(items)
  x=items(k); off=x.offset;
  if strcmp(x.domain,'IDQ'), off=120; end
  % Identical paths are deduplicated here. Equivalent whole-session IDQ
  % representations are checked after extraction, allowing a valid legacy
  % source to be used if the preferred canonical source cannot be read.
  key=x.source;
  if isKey(seen,key)
    keep(k)=false; notes(end+1,:)={x.domain,x.source,'Duplicate source representation; preferred canonical source retained'}; %#ok<AGROW>
  else, seen(key)=true;
  end
end
items=items(keep);
end
