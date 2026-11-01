function [cache,manifest] = prepareLocalCache(files,folder,cfg,rebuild)
%PREPARELOCALCACHE Open only new/changed large source files, one at a time.
if nargin<4, rebuild=false; end
if ~isfolder(folder), mkdir(folder); end
files=cellstr(string(files(:))); parts=cell(numel(files),1); rows=cell(numel(files),6);
extraction=struct('version',"IDS-window-cache-v2",'animal',cfg.animal, ...
  'preWindowMS',cfg.preWindowMS,'stepWindowMS',cfg.stepWindowMS, ...
  'keepRawStepNoise',cfg.keepRawStepNoise);
for k=1:numel(files)
  d=dir(files{k}); assert(isscalar(d),'Missing or ambiguous source: %s',files{k});
  path=fullfile(d.folder,d.name);
  fingerprint=struct('path',string(path),'bytes',d.bytes,'datenum',d.datenum,'extraction',extraction);
  cacheFile=fullfile(folder,[pathKey(path) '.mat']); reuse=false;
  if ~rebuild && isfile(cacheFile)
    S=load(cacheFile,'fingerprint');
    reuse=isequaln(S.fingerprint,fingerprint);
  end
  if reuse
    S=load(cacheFile,'fragment'); fragment=S.fragment; status="reused";
  else
    fprintf('Extracting %d/%d: %s (%.1f MB)\n',k,numel(files),d.name,d.bytes/2^20);
    fragment=mtpnorm.readIDS({path},cfg);
    after=dir(path);
    assert(after.bytes==d.bytes && after.datenum==d.datenum, ...
      'Source changed during extraction; retry: %s',path);
    tmp=[cacheFile '.tmp.mat']; save(tmp,'fragment','fingerprint','-v7.3');
    movefile(tmp,cacheFile,'f'); status="extracted";
  end
  parts{k}=fragment;
  rows(k,:)={string(d.name),d.bytes,d.datenum,status,height(fragment.trials),string(cacheFile)};
end
assert(~isempty(parts),'No source files selected.');
cache=parts{1}; names={'trials','preCohPC','stepCohPC','streamDirDeg','sourceAudit'};
for k=1:numel(names)
  values=cellfun(@(s)s.(names{k}),parts,'UniformOutput',false);
  cache.(names{k})=vertcat(values{:});
end
if cfg.keepRawStepNoise
  values=cellfun(@(s)s.rawStepNoisePC,parts,'UniformOutput',false);
  cache.rawStepNoisePC=vertcat(values{:});
else
  cache=rmfield(cache,'rawStepNoisePC');
end
cache.config=cfg;
keys=string(cache.trials.sessionID)+"/"+string(cache.trials.parentTrialIdx);
assert(numel(unique(keys))==numel(keys), ...
  'Duplicate parent trials: select either parent sessions or ProbeSessions, not both.');
manifest=cell2table(rows,'VariableNames',{'fileName','sourceBytes','sourceDatenum', ...
  'cacheStatus','nValidTrials','fragmentPath'});
fprintf('Cache: %d reused, %d extracted; %d trials.\n', ...
  nnz(manifest.cacheStatus=="reused"),nnz(manifest.cacheStatus=="extracted"),height(cache.trials));
end

function key=pathKey(path)
% Portable filename identifier; the full fingerprint is always verified.
bytes=double(unicode2native(char(path),'UTF-8')); a=5381; b=2166136261;
for k=1:numel(bytes)
  a=mod(a*33+bytes(k),4294967291);
  b=mod(b*65599+bytes(k),4294967279);
end
key=sprintf('%08x%08x',uint32(a),uint32(b));
end
