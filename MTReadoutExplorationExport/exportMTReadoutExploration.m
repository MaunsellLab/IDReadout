function bundlePath = exportMTReadoutExploration(varargin)
%EXPORTMTREADOUTEXPLORATION Export compact inputs for automated model search.
% exportMTReadoutExploration  -- choose the common IDReadout project root once.
% Does not fit a model, optimize noise, or change existing analysis products.
p=inputParser;
p.addParameter('ProjectRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('OutputFolder','',@(x)ischar(x)||isstring(x));
p.addParameter('Animal','Meetz',@(x)ischar(x)||isstring(x));
p.addParameter('Domains',{'IDS','IDQ','IDJ','IDR'},@(x)iscellstr(x)||isstring(x));
p.addParameter('BinMS',50,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('ChooseFolder',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('Rebuild',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('MaxFiles',Inf,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&(isinf(x)||x==fix(x)));
p.addParameter('IncludeAnalyses',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); O=p.Results;
O.Domains=cellstr(upper(string(O.Domains)));
root=char(O.ProjectRoot);
if isempty(root)&&~O.ChooseFolder&&ispref('MTReadoutExploration','projectRoot')
  root=getpref('MTReadoutExploration','projectRoot');
end
if O.ChooseFolder||isempty(root)||~isfolder(root)
  start=pwd;
  if ispref('MTReadoutPNorm','dataRoot')
    start=getpref('MTReadoutPNorm','dataRoot');
    for k=1:3
      if isfolder(fullfile(start,'IDS'))||isfolder(fullfile(start,'IDQ')), break; end
      parent=fileparts(start); if isempty(parent), break; end; start=parent;
    end
  end
  chosen=uigetdir(start,'Choose the common IDReadout root containing IDS, IDQ and IDJ');
  assert(~isequal(chosen,0),'Folder selection cancelled.'); root=chosen;
end
out=char(O.OutputFolder); if isempty(out), out=fullfile(root,'MTReadoutExplorationOutput'); end
if ~isfolder(out), mkdir(out); end
for folder={'sources','analyses','analysis_code'}
  if ~isfolder(fullfile(out,folder{1})), mkdir(fullfile(out,folder{1})); end
end
[items,notes]=mtexport.discover(root,O);
assert(~isempty(items),'No source files found. Choose the common project root and check the report/exclusions.');
setpref('MTReadoutExploration','projectRoot',root);
subset=numel(items)>O.MaxFiles;
if subset, items=items(1:O.MaxFiles); end
files={}; rows=cell(0,11); keys=containers.Map('KeyType','char','ValueType','logical');
for k=1:numel(items)
  item=items(k); relative=fullfile('sources',['source_' mtexport.id(item.source) '.mat']);
  path=fullfile(out,relative); d=dir(item.source); c=[];
  if ~isempty(item.compact), c=dir(item.compact); end
  fp=struct('version',3,'source',item.source,'bytes',d.bytes,'datenum',d.datenum, ...
    'compact',item.compact,'compactBytes',[],'compactDatenum',[],'binMS',O.BinMS, ...
    'domain',item.domain,'animal',item.animal,'randomBase',item.randomBase);
  if ~isempty(c), fp.compactBytes=c.bytes; fp.compactDatenum=c.datenum; end
  try
    reuse=false;
    if ~O.Rebuild&&isfile(path)
      A=load(path,'fingerprint'); reuse=isfield(A,'fingerprint')&&isequaln(A.fingerprint,fp);
    end
    if reuse
      A=load(path,'source'); source=A.source; state='reused';
    else
      fprintf('Exporting %d/%d: %s (%s)\n',k,numel(items),item.domain,item.source);
      source=mtexport.readSource(item,O.BinMS); %#ok<NASGU>
      after=dir(item.source);
      assert(after.bytes==d.bytes&&after.datenum==d.datenum,'Source changed during export; retry.');
      fingerprint=fp; %#ok<NASGU>
      tmp=[path '.tmp.mat']; save(tmp,'source','fingerprint','-v7'); movefile(tmp,path,'f'); state='exported';
    end
    % Check identities before including a shard, rather than double-counting
    % duplicated parent/derived trials in a later behavioral likelihood.
    existing=false(numel(source.parent_trial_index),1);
    for j=1:numel(source.parent_trial_index)
      key=sprintf('%s|%s|%s|%d',source.meta.domain,source.meta.animal, ...
        source.meta.session_id,source.parent_trial_index(j));
      existing(j)=isKey(keys,key);
    end
    if all(existing)
      notes(end+1,:)={item.domain,item.source,'Duplicate trial representation; earlier preferred source retained'}; %#ok<AGROW>
      continue
    end
    assert(~any(existing),'Partially overlapping source representations for %s.',source.meta.session_id);
    for j=1:numel(source.parent_trial_index)
      key=sprintf('%s|%s|%s|%d',source.meta.domain,source.meta.animal, ...
        source.meta.session_id,source.parent_trial_index(j)); keys(key)=true;
    end
    files{end+1,1}=relative; %#ok<AGROW>
    rows(end+1,:)={source.meta.domain,source.meta.animal,source.meta.session_id, ...
      source.meta.source_name,strrep(relative,filesep,'/'),source.meta.n_trials, ...
      source.meta.n_fit_trials,source.meta.n_excluded_topup, ...
      mat2str(unique(source.probe_offset_deg(isfinite(source.probe_offset_deg)))'), ...
      source.meta.input_representation,state}; %#ok<AGROW>
  catch ME
    if strcmp(ME.identifier,'MTReadoutExport:NoEligibleTrials')
      notes(end+1,:)={item.domain,item.source,['Excluded: ' ME.message]}; %#ok<AGROW>
    else
      notes(end+1,:)={item.domain,item.source,['Export error: ' ME.message]}; %#ok<AGROW>
      warning('MTReadoutExploration:Source','%s',ME.message);
    end
  end
end
assert(~isempty(rows),'No sources exported; inspect source metadata and errors above.');
manifest=cell2table(rows,'VariableNames',{'domain','animal','session_id','source_name', ...
  'bundle_file','n_trials','n_fit_trials','n_excluded_topup','probe_offsets_deg', ...
  'input_representation','cache_status'});
writetable(manifest,fullfile(out,'SourceManifest.csv'));
report=cell2table(notes,'VariableNames',{'domain','source','note'});
writetable(report,fullfile(out,'ExportNotes.csv'));
if O.IncludeAnalyses
  [extra,extraNotes]=mtexport.exportAnalyses(root,out,O); files=[files;extra(:)];
  report=[report;cell2table(extraNotes,'VariableNames',report.Properties.VariableNames)];
  writetable(report,fullfile(out,'ExportNotes.csv'));
end
coverage=cell(numel(O.Domains),4);
for k=1:numel(O.Domains)
  mask=strcmpi(manifest.domain,O.Domains{k});
  coverage(k,:)={O.Domains{k},nnz(mask),sum(manifest.n_trials(mask)),sum(manifest.n_fit_trials(mask))};
end
coverage=cell2table(coverage,'VariableNames',{'domain','n_files','n_trials','n_fit_trials'});
writetable(coverage,fullfile(out,'Coverage.csv'));
schema=struct('schema_version',1,'animal_selection',char(O.Animal), ...
  'requested_domains',{O.Domains},'bin_ms',O.BinMS,'subset_run',subset, ...
  'created_at',char(datetime('now')),'pooling_exponent','not imposed by exporter', ...
  'axis_order','trial, physical RDK, stream, time bin','side_index_base',1, ...
  'noise_units','coherence percentage points','idq_selection','randomBaseDir=true only', ...
  'topup_selection','excluded','n_export_errors',nnz(contains(string(report.note),'error:')));
fid=fopen(fullfile(out,'Schema.json'),'w'); assert(fid>=0); fprintf(fid,'%s',jsonencode(schema)); fclose(fid);
fid=fopen(fullfile(out,'README.txt'),'w'); assert(fid>=0);
fprintf(fid,'MT readout exploration INPUT bundle, schema 1\n');
fprintf(fid,'Plain numeric MAT v7 shards, CSV manifest/coverage, saved analysis metadata and code.\n');
fprintf(fid,'No population model is fitted and no p is imposed. Only randomized-base IDQ is exported.\n');
fprintf(fid,'Check Coverage.csv and ExportNotes.csv: missing domains/errors mean partial coverage.\n');
fprintf(fid,'noise_pc contains signed fluctuations; add deterministic_pre_pc or deterministic_step_pc by epoch.\n');
fprintf(fid,'Inactive padded streams have stream_mask=false and NaN inputs. Side indices are 1 or 2.\n');
fprintf(fid,'Window totals are stored directly, in addition to binned full-trial noise.\n');
fprintf(fid,'Raw allocation fractions are preserved; near-equal thirds also receive nominal equal-third labels.\n');
fprintf(fid,'Current source shards are cached locally; the ZIP explicitly includes only this selection.\n'); fclose(fid);
files=[files;{'SourceManifest.csv';'Coverage.csv';'ExportNotes.csv';'Schema.json';'README.txt'}];
bundlePath=fullfile(out,'MTReadoutExplorationInputs.zip'); zip(bundlePath,files,out);
disp(coverage);
fprintf('\nExported %d sources. Send: %s\n',height(manifest),bundlePath);
if any(coverage.n_files==0)||schema.n_export_errors>0
  fprintf('Coverage is partial; inspect Coverage.csv / ExportNotes.csv. The bundle records these limitations.\n');
end
end
