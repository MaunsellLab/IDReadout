function [cache,report] = syncCache(scanManifest,readerFcn,cachePath,varargin)
%SYNCCACHE Incrementally synchronize a manifest-backed trial cache.
%
% Unchanged source MAT files are never opened. New files are appended;
% modified files are re-extracted and replace only their prior rows. Missing
% files stop the update by default. Each extracted ProbeSession is saved as
% an atomic restart checkpoint before the next source is opened. The final
% consolidated cache is also saved atomically.

p=inputParser;
p.addParameter('ExtractorVersion',"ProbeSessionAdapter-v4-full-prestep",@(x)ischar(x)||isstring(x));
p.addParameter('Rebuild',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('DryRun',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('MissingFilePolicy',"error",@(x)ismember(string(x),["error" "retain" "remove"]));
p.addParameter('RecomputeKernels',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('ClearCheckpoints',false,@(x)islogical(x)&&isscalar(x));
p.addParameter('CheckpointFolder',"",@(x)ischar(x)||isstring(x));
p.addParameter('Verbose',true,@(x)islogical(x)&&isscalar(x));
p.parse(varargin{:}); opts=p.Results;
extractorVersion=string(opts.ExtractorVersion);
checkpointFolder=string(opts.CheckpointFolder);
if strlength(checkpointFolder)==0, checkpointFolder=string(cachePath)+".staging"; end
validateScanManifest(scanManifest);
scanManifest.sourceFile=string(scanManifest.sourceFile);

hasExisting=isfile(cachePath)&&~opts.Rebuild;
if hasExisting
  S=load(cachePath,'cache'); assert(isfield(S,'cache'),'Cache file lacks cache variable.');
  old=S.cache;
  assert(isfield(old,'schemaVersion') && ...
    string(old.schemaVersion)=="MaxPoolingTrialCache-v2", ...
    ['Existing cache predates full prestep-stream adaptation. Rerun with ' ...
     'Rebuild=true; the previous final cache will remain untouched until ' ...
     'the complete replacement has been validated and saved.']);
  maxpoolfit.validateCache(old);
  assert(isfield(old,'sourceManifest'), ...
    ['Existing cache predates incremental manifests. Run once with ' ...
     'Rebuild=true to create the incremental cache.']);
  oldManifest=old.sourceManifest;
else
  old=[]; oldManifest=emptyManifest();
end

[newPaths,changedPaths,missingPaths,unchangedPaths]=compareManifests( ...
  scanManifest,oldManifest,extractorVersion);
if opts.Rebuild
  newPaths=scanManifest.sourceFile; changedPaths=strings(0,1);
  missingPaths=strings(0,1); unchangedPaths=strings(0,1);
end
refreshPaths=[newPaths;changedPaths];
checkpointIndex=emptyCheckpointIndex();
if ~opts.ClearCheckpoints
  checkpointIndex=loadCheckpointIndex(checkpointFolder,opts.Verbose);
end
checkpointPaths=findReusableCheckpointPaths( ...
  refreshPaths,scanManifest,checkpointIndex,extractorVersion);
report=makeReport(newPaths,changedPaths,missingPaths,unchangedPaths, ...
  checkpointPaths,opts,checkpointFolder);
if opts.Verbose, printReport(report); end
if opts.DryRun
  if hasExisting, cache=old; else, cache=[]; end
  if hasExisting, report.nCachedTrials=height(old.trials); else, report.nCachedTrials=0; end
  return
end
if ~isempty(missingPaths)&&string(opts.MissingFilePolicy)=="error"
  error('maxpoolfit:syncCache:MissingSources', ...
    ['%d cached source files are absent from the scan. No changes were made. ' ...
     'Use MissingFilePolicy="retain" or "remove" after reviewing report.missingFiles.'], ...
    numel(missingPaths));
end
if opts.ClearCheckpoints
  clearCheckpointFiles(checkpointFolder,opts.Verbose);
  checkpointIndex=emptyCheckpointIndex();
  checkpointPaths=strings(0,1);
end
if hasExisting&&isempty(newPaths)&&isempty(changedPaths)&&isempty(missingPaths)
  cache=old; report.nCachedTrials=height(cache.trials); return
end

if opts.Rebuild
  base=[];
else
  removePaths=changedPaths;
  if string(opts.MissingFilePolicy)=="remove", removePaths=[removePaths;missingPaths]; end
  base=removeSourceFiles(old,removePaths);
end

partFiles=strings(numel(refreshPaths),1); nCheckpointUsed=0;
for k=1:numel(refreshPaths)
  scanRow=scanManifest(scanManifest.sourceFile==refreshPaths(k),:);
  assert(height(scanRow)==1,'Directory manifest did not resolve one row for %s.',refreshPaths(k));
  verifyFileFingerprint(refreshPaths(k),scanRow);
  checkpointFile=findReusableCheckpoint( ...
    refreshPaths(k),scanRow,checkpointIndex,extractorVersion);
  fragment=[];
  if strlength(checkpointFile)>0
    try
      fragment=loadCheckpointFragment(checkpointFile,refreshPaths(k),extractorVersion);
      validateFragment(fragment,refreshPaths(k),extractorVersion);
      nCheckpointUsed=nCheckpointUsed+1;
      if opts.Verbose
        fprintf('Resuming %d/%d from checkpoint: %s\n', ...
          k,numel(refreshPaths),refreshPaths(k));
      end
    catch ME
      warning('maxpoolfit:syncCache:InvalidCheckpoint', ...
        'Ignoring invalid checkpoint for %s: %s',refreshPaths(k),ME.message);
      fragment=[]; checkpointFile="";
    end
  end
  if isempty(fragment)
    if opts.Verbose
      fprintf('Extracting %d/%d: %s\n',k,numel(refreshPaths),refreshPaths(k));
    end
    fragment=readerFcn(char(refreshPaths(k)));
    verifyFileFingerprint(refreshPaths(k),scanRow);
    validateFragment(fragment,refreshPaths(k),extractorVersion);
    [checkpointFile,checkpointRow]=saveFragmentCheckpoint( ...
      checkpointFolder,fragment,refreshPaths(k),scanRow,extractorVersion);
    checkpointIndex=[checkpointIndex;checkpointRow]; %#ok<AGROW>
    if opts.Verbose, fprintf('  checkpoint saved\n'); end
  end
  partFiles(k)=checkpointFile;
  clear fragment
end
cache=combineBaseAndPartFiles(base,partFiles);

currentManifest=scanManifest;
currentManifest.extractorVersion=repmat(extractorVersion,height(currentManifest),1);
if string(opts.MissingFilePolicy)=="retain"&&~isempty(missingPaths)
  keepOld=ismember(string(oldManifest.sourceFile),missingPaths);
  retained=oldManifest(keepOld,intersect(currentManifest.Properties.VariableNames, ...
    oldManifest.Properties.VariableNames,'stable'));
  currentManifest=[currentManifest;retained];
end
currentManifest.nCachedTrials=zeros(height(currentManifest),1);
for k=1:height(currentManifest)
  currentManifest.nCachedTrials(k)=sum(string(cache.trials.sourceFile)==currentManifest.sourceFile(k));
end
cache.sourceManifest=sortrows(currentManifest,'sourceFile');
cache.sourceFiles=cache.sourceManifest.sourceFile;
cache.schemaVersion="MaxPoolingTrialCache-v2";
cache.incrementalCacheVersion="IncrementalCache-v1";
cache.extractorVersion=extractorVersion;
cache.updatedAt=datetime('now');
if ~isfield(cache,'createdAt'), cache.createdAt=cache.updatedAt; end
cache.createdBy=mfilename;
historyEntry=struct('updatedAt',cache.updatedAt,'nNew',numel(newPaths), ...
  'nChanged',numel(changedPaths),'nMissing',numel(missingPaths), ...
  'nUnchanged',numel(unchangedPaths),'extractorVersion',extractorVersion);
if isfield(cache,'updateHistory'), cache.updateHistory(end+1)=historyEntry;
else, cache.updateHistory=historyEntry; end
if opts.RecomputeKernels
  [cache,cache.kernels]=maxpoolfit.crossfitKernels(cache);
else
  if isfield(cache,'kernels'), cache=rmfield(cache,'kernels'); end
  cache.kernelStatus="stale";
end
maxpoolfit.validateCache(cache);
atomicSave(cachePath,cache);
report.nCachedTrials=height(cache.trials);
report.nCheckpointUsed=nCheckpointUsed;
try
  clearCheckpointFiles(checkpointFolder,opts.Verbose);
catch ME
  warning('maxpoolfit:syncCache:CheckpointCleanupFailed', ...
    'Final cache was saved, but checkpoint cleanup failed: %s',ME.message);
end
end

function validateScanManifest(M)
required={'sourceFile','domain','animal','bytes','modifiedDatenum'};
assert(istable(M)&&all(ismember(required,M.Properties.VariableNames)), ...
  'scanManifest lacks required variables.');
assert(numel(unique(string(M.sourceFile)))==height(M),'scanManifest paths must be unique.');
end

function verifyFileFingerprint(path,row)
D=dir(char(path));
assert(isscalar(D)&&~D.isdir,'Source disappeared or became ambiguous during update: %s',path);
assert(double(D.bytes)==double(row.bytes)&&double(D.datenum)==double(row.modifiedDatenum), ...
  ['Source changed after the directory scan: %s. No cache was written; ' ...
   'rerun the synchronization.'],path);
end

function M=emptyManifest()
M=table(strings(0,1),strings(0,1),strings(0,1),zeros(0,1),zeros(0,1), ...
  strings(0,1),zeros(0,1),'VariableNames',{'sourceFile','domain','animal', ...
  'bytes','modifiedDatenum','extractorVersion','nCachedTrials'});
end

function [newP,changedP,missingP,unchangedP]=compareManifests(scan,old,version)
scanP=string(scan.sourceFile); oldP=string(old.sourceFile);
newP=setdiff(scanP,oldP,'stable'); missingP=setdiff(oldP,scanP,'stable');
common=intersect(scanP,oldP,'stable'); changedP=strings(0,1); unchangedP=strings(0,1);
for k=1:numel(common)
  i=find(scanP==common(k),1); j=find(oldP==common(k),1);
  changed=double(scan.bytes(i))~=double(old.bytes(j)) || ...
    double(scan.modifiedDatenum(i))~=double(old.modifiedDatenum(j)) || ...
    string(scan.domain(i))~=string(old.domain(j)) || ...
    string(scan.animal(i))~=string(old.animal(j)) || ...
    ~ismember('extractorVersion',old.Properties.VariableNames) || ...
    string(old.extractorVersion(j))~=version;
  if changed, changedP(end+1,1)=common(k); %#ok<AGROW>
  else, unchangedP(end+1,1)=common(k); %#ok<AGROW>
  end
end
end

function cache=removeSourceFiles(cache,paths)
if isempty(paths), return; end
keep=~ismember(string(cache.trials.sourceFile),string(paths));
cache.trials=cache.trials(keep,:); cache.streamDirDeg=cache.streamDirDeg(keep,:);
cache.stepCohPC=cache.stepCohPC(keep,:); cache.driftDirDeg=cache.driftDirDeg(keep,:);
cache.preCohPC=cache.preCohPC(keep,:);
cache.preDriftCohPC=cache.preDriftCohPC(keep,:);
if isfield(cache,'raw')
  cache.raw.stepStreamCohPC=cache.raw.stepStreamCohPC(keep,:);
  cache.raw.kernelDriverByFrame=cache.raw.kernelDriverByFrame(keep,:);
end
if isfield(cache,'extractionAudits')
  auditPaths=strings(numel(cache.extractionAudits),1);
  for k=1:numel(cache.extractionAudits)
    auditPaths(k)=string(cache.extractionAudits{k}.sourcePath);
  end
  cache.extractionAudits=cache.extractionAudits(~ismember(auditPaths,string(paths)));
end
end

function cache=combineBaseAndPartFiles(base,partFiles)
if isempty(base)
  assert(~isempty(partFiles),'No source files available to create cache.');
  S=load(char(partFiles(1)),'fragment');
  assert(isfield(S,'fragment'),'Checkpoint lacks fragment: %s',partFiles(1));
  cache=S.fragment; start=2;
  audits={}; if isfield(cache,'extractionAudit'), audits={cache.extractionAudit}; end
else
  cache=base; start=1;
  if isfield(cache,'extractionAudits'), audits=cache.extractionAudits; else, audits={}; end
end
for k=start:numel(partFiles)
  S=load(char(partFiles(k)),'fragment');
  assert(isfield(S,'fragment'),'Checkpoint lacks fragment: %s',partFiles(k));
  P=S.fragment;
  cache.trials=[cache.trials;P.trials];
  cache.streamDirDeg=[cache.streamDirDeg;P.streamDirDeg];
  cache.stepCohPC=[cache.stepCohPC;P.stepCohPC];
  cache.preCohPC=[cache.preCohPC;P.preCohPC];
  cache.driftDirDeg=[cache.driftDirDeg;P.driftDirDeg];
  cache.preDriftCohPC=[cache.preDriftCohPC;P.preDriftCohPC];
  cache.raw.stepStreamCohPC=[cache.raw.stepStreamCohPC;P.raw.stepStreamCohPC];
  cache.raw.kernelDriverByFrame=[cache.raw.kernelDriverByFrame;P.raw.kernelDriverByFrame];
  if isfield(P,'extractionAudit'), audits{end+1,1}=P.extractionAudit; end %#ok<AGROW>
end
if isempty(base)&&isfield(cache,'extractionAudit'), cache=rmfield(cache,'extractionAudit'); end
cache.extractionAudits=audits;
end

function validateFragment(fragment,sourcePath,extractorVersion)
assert(isfield(fragment,'extractorVersion')&& ...
  string(fragment.extractorVersion)==extractorVersion, ...
  'Reader extractor version does not match requested version.');
maxpoolfit.validateCache(partialCache(fragment),'AllowEmpty',true);
assert(ismember('sourceFile',fragment.trials.Properties.VariableNames)&& ...
  all(string(fragment.trials.sourceFile)==sourcePath), ...
  'Reader returned trials with the wrong sourceFile for %s.',sourcePath);
end

function I=emptyCheckpointIndex()
I=table(strings(0,1),zeros(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
  'VariableNames',{'sourceFile','bytes','modifiedDatenum', ...
  'extractorVersion','checkpointFile'});
end

function I=loadCheckpointIndex(folder,verbose)
I=emptyCheckpointIndex();
if ~isfolder(folder), return; end
D=dir(char(fullfile(folder,'checkpoint_*.mat')));
for k=1:numel(D)
  checkpointFile=string(fullfile(D(k).folder,D(k).name));
  try
    S=load(char(checkpointFile),'checkpointMetadata');
    assert(isfield(S,'checkpointMetadata'),'missing checkpointMetadata');
    M=S.checkpointMetadata;
    required={'formatVersion','sourceFile','bytes','modifiedDatenum','extractorVersion'};
    assert(isstruct(M)&&all(isfield(M,required)),'incomplete checkpointMetadata');
    assert(string(M.formatVersion)=="MaxPoolingFragmentCheckpoint-v1", ...
      'unsupported checkpoint format');
    row=table(string(M.sourceFile),double(M.bytes),double(M.modifiedDatenum), ...
      string(M.extractorVersion),checkpointFile,'VariableNames', ...
      I.Properties.VariableNames);
    I=[I;row]; %#ok<AGROW>
  catch ME
    if verbose
      warning('maxpoolfit:syncCache:UnreadableCheckpoint', ...
        'Ignoring checkpoint %s: %s',checkpointFile,ME.message);
    end
  end
end
end

function paths=findReusableCheckpointPaths(refreshPaths,scanManifest,I,version)
paths=strings(0,1);
for k=1:numel(refreshPaths)
  row=scanManifest(scanManifest.sourceFile==refreshPaths(k),:);
  if strlength(findReusableCheckpoint(refreshPaths(k),row,I,version))>0
    paths(end+1,1)=refreshPaths(k); %#ok<AGROW>
  end
end
end

function checkpointFile=findReusableCheckpoint(sourcePath,scanRow,I,version)
checkpointFile="";
if isempty(I), return; end
match=string(I.sourceFile)==sourcePath & ...
  double(I.bytes)==double(scanRow.bytes) & ...
  double(I.modifiedDatenum)==double(scanRow.modifiedDatenum) & ...
  string(I.extractorVersion)==version;
idx=find(match,1,'last');
if ~isempty(idx), checkpointFile=I.checkpointFile(idx); end
end

function fragment=loadCheckpointFragment(checkpointFile,sourcePath,version)
S=load(char(checkpointFile),'checkpointMetadata','fragment');
assert(isfield(S,'checkpointMetadata')&&isfield(S,'fragment'), ...
  'Checkpoint is incomplete: %s',checkpointFile);
assert(string(S.checkpointMetadata.sourceFile)==sourcePath&& ...
  string(S.checkpointMetadata.extractorVersion)==version, ...
  'Checkpoint identity does not match its requested source.');
fragment=S.fragment;
end

function [checkpointFile,row]=saveFragmentCheckpoint( ...
    folder,fragment,sourcePath,scanRow,version)
if ~isfolder(folder), mkdir(char(folder)); end
[~,token]=fileparts(tempname(char(folder)));
checkpointFile=string(fullfile(folder,"checkpoint_"+string(token)+".mat"));
partialFile=char(fullfile(folder,"partial_"+string(token)+".mat"));
cleanup=onCleanup(@()deleteIfPresent(partialFile)); %#ok<NASGU>
checkpointMetadata=struct('formatVersion',"MaxPoolingFragmentCheckpoint-v1", ...
  'sourceFile',sourcePath,'bytes',double(scanRow.bytes), ...
  'modifiedDatenum',double(scanRow.modifiedDatenum), ...
  'extractorVersion',version,'createdAt',datetime('now'));
save(partialFile,'checkpointMetadata','fragment','-v7.3');
[ok,msg]=movefile(partialFile,char(checkpointFile),'f');
assert(ok,'Atomic checkpoint save failed: %s',msg);
row=table(sourcePath,double(scanRow.bytes),double(scanRow.modifiedDatenum), ...
  version,checkpointFile,'VariableNames',checkpointVariableNames());
end

function names=checkpointVariableNames()
names={'sourceFile','bytes','modifiedDatenum','extractorVersion','checkpointFile'};
end

function clearCheckpointFiles(folder,verbose)
if ~isfolder(folder), return; end
D=[dir(char(fullfile(folder,'checkpoint_*.mat'))); ...
   dir(char(fullfile(folder,'partial_*.mat')))];
for k=1:numel(D)
  delete(fullfile(D(k).folder,D(k).name));
end
remaining=dir(char(folder));
remaining=remaining(~ismember({remaining.name},{'.','..'}));
if isempty(remaining), rmdir(char(folder)); end
if verbose, fprintf('Restart checkpoints cleared after successful cache save.\n'); end
end

function C=partialCache(P)
C=P; C.schemaVersion="MaxPoolingTrialCache-v2";
end

function R=makeReport(newP,changedP,missingP,unchangedP,checkpointP,opts,checkpointFolder)
R=struct('newFiles',newP,'changedFiles',changedP,'missingFiles',missingP, ...
  'unchangedFiles',unchangedP,'nNew',numel(newP),'nChanged',numel(changedP), ...
  'nMissing',numel(missingP),'nUnchanged',numel(unchangedP), ...
  'checkpointFiles',checkpointP,'nCheckpoint',numel(checkpointP), ...
  'checkpointFolder',checkpointFolder,'dryRun',opts.DryRun, ...
  'missingFilePolicy',string(opts.MissingFilePolicy));
files=[newP;changedP;missingP;unchangedP];
actions=[repmat("extract",numel(newP),1);repmat("replace",numel(changedP),1); ...
  repmat("missing",numel(missingP),1);repmat("reuse",numel(unchangedP),1)];
actions(ismember(files,checkpointP))="resume";
R.actions=table(files,actions,'VariableNames',{'sourceFile','action'});
end

function printReport(R)
fprintf('Cache scan: %d new, %d changed, %d missing, %d unchanged.\n', ...
  R.nNew,R.nChanged,R.nMissing,R.nUnchanged);
if R.nCheckpoint>0
  fprintf('Restart checkpoints available for %d source files.\n',R.nCheckpoint);
end
end

function atomicSave(path,cache)
folder=fileparts(path); if isempty(folder), folder=pwd; end
if ~isfolder(folder), mkdir(folder); end
tmp=[tempname(folder) '.mat']; cleanup=onCleanup(@()deleteIfPresent(tmp));
save(tmp,'cache','-v7.3');
[ok,msg]=movefile(tmp,path,'f'); assert(ok,'Atomic cache replacement failed: %s',msg);
clear cleanup
end

function deleteIfPresent(path)
if isfile(path), delete(path); end
end
