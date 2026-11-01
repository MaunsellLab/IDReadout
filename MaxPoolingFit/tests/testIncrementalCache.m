function tests=testIncrementalCache
tests=functiontests(localfunctions);
end

function testDiscoveryUsesCentralHardExcludes(testCase)
root=tempname; cleanup=onCleanup(@()removeTree(root)); %#ok<NASGU>
projectCodeRoot=fileparts(fileparts(fileparts(mfilename('fullpath'))));
oldPath=path; addpath(projectCodeRoot); pathCleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
mkdir(fullfile(root,'IDR','Data','Probe10','ProbeSessions'));
folder=fullfile(root,'IDR','Data','Probe25','ProbeSessions'); mkdir(folder);
excluded=fullfile(folder,'IDReadout_Meetz_20260114_Probe25.mat'); touch(excluded,'a');
included=fullfile(folder,'IDReadout2_Meetz_20260615_Probe25.mat'); touch(included,'a');
M=maxpoolfit.discoverProbeSessions(root,'Domains',"IDR",'Animals',"Meetz");
verifyEqual(testCase,height(M),1);
verifyEqual(testCase,M.sourceFile,string(included));
end

function testInterruptedBuildResumesFromPerSessionCheckpoint(testCase)
root=tempname; mkdir(root); cleanup=onCleanup(@()removeTree(root)); %#ok<NASGU>
file1=fullfile(root,'IDReadout2_Meetz_A_Probe25.mat'); touch(file1,'a');
file2=fullfile(root,'IDReadout2_Meetz_B_Probe25.mat'); touch(file2,'a');
M=makeManifest([string(file1);string(file2)]);
cachePath=fullfile(root,'cache.mat'); checkpointFolder=[cachePath '.staging'];
readerCalls=strings(0,1); failureEnabled=true;

verifyError(testCase,@()maxpoolfit.syncCache(M,@checkpointReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false),'test:InjectedReaderFailure');
verifyFalse(testCase,isfile(cachePath));
D=dir(fullfile(checkpointFolder,'checkpoint_*.mat'));
verifyEqual(testCase,numel(D),1);
[~,preview]=maxpoolfit.syncCache(M,@checkpointReader,cachePath, ...
  'DryRun',true,'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,preview.nCheckpoint,1);
verifyEqual(testCase,sum(preview.actions.action=="resume"),1);

failureEnabled=false;
[C,R]=maxpoolfit.syncCache(M,@checkpointReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,height(C.trials),2);
verifyEqual(testCase,R.nCheckpointUsed,1);
verifyEqual(testCase,sum(readerCalls==string(file1)),1);
verifyEqual(testCase,sum(readerCalls==string(file2)),2);
verifyFalse(testCase,isfolder(checkpointFolder));

  function F=checkpointReader(path)
    readerCalls(end+1,1)=string(path); %#ok<AGROW>
    if failureEnabled && string(path)==string(file2)
      error('test:InjectedReaderFailure','Injected failure after one checkpoint.');
    end
    F=mockReader(path);
  end
end

function testNewChangedUnchangedAndDryRunMissing(testCase)
root=tempname; mkdir(root); cleanup=onCleanup(@()removeTree(root)); %#ok<NASGU>
file1=fullfile(root,'IDReadout2_Meetz_A_Probe25.mat'); touch(file1,'a');
cachePath=fullfile(root,'cache.mat');

M=makeManifest(file1);
[C,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,R.nNew,1); verifyEqual(testCase,height(C.trials),1);

[C,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,R.nUnchanged,1); verifyEqual(testCase,height(C.trials),1);

touch(file1,'ab'); M=makeManifest(file1);
[C,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,R.nChanged,1); verifyEqual(testCase,height(C.trials),1);

file2=fullfile(root,'IDReadout2_Meetz_B_Probe25.mat'); touch(file2,'a');
M=makeManifest([string(file1);string(file2)]);
[C,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,R.nNew,1); verifyEqual(testCase,height(C.trials),2);

delete(file2); M=makeManifest(file1);
[~,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'DryRun',true,'Verbose',false);
verifyEqual(testCase,R.nMissing,1);
verifyEqual(testCase,R.actions.action(R.actions.action=="missing"),"missing");
verifyError(testCase,@()maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'RecomputeKernels',false,'Verbose',false), ...
  'maxpoolfit:syncCache:MissingSources');
[C,R]=maxpoolfit.syncCache(M,@mockReader,cachePath, ...
  'MissingFilePolicy','remove','RecomputeKernels',false,'Verbose',false);
verifyEqual(testCase,R.nMissing,1); verifyEqual(testCase,height(C.trials),1);
end

function F=mockReader(path)
F=maxpoolfit.syntheticCache('TrialsPerCondition',1);
F.trials=F.trials(1,:); F.streamDirDeg=F.streamDirDeg(1,:);
F.stepCohPC=F.stepCohPC(1,:); F.preCohPC=F.preCohPC(1,:);
F.driftDirDeg=F.driftDirDeg(1,:);
F.preDriftCohPC=F.preDriftCohPC(1,:);
[~,name]=fileparts(path); F.trials.sessionID=string(name);
F.trials.sourceFile=string(path);
F.raw.stepStreamCohPC={repmat(F.stepCohPC{1,1}',2,1), ...
  repmat(F.stepCohPC{1,2}',2,1)};
F.raw.kernelDriverByFrame={[1;-1]};
F.extractorVersion="ProbeSessionAdapter-v4-full-prestep";
F.extractionAudit=struct('sourcePath',string(path),'domain',"IDR", ...
  'extractorVersion',F.extractorVersion);
end

function M=makeManifest(paths)
paths=string(paths); paths=paths(:);
n=numel(paths); bytes=zeros(n,1); modified=zeros(n,1);
for k=1:n
  D=dir(char(paths(k)));
  assert(isscalar(D)&&~D.isdir,'Test source file was not found: %s',paths(k));
  bytes(k)=D.bytes; modified(k)=D.datenum;
end
M=table(paths,repmat("IDR",n,1),repmat("Meetz",n,1),bytes,modified, ...
  'VariableNames',{'sourceFile','domain','animal','bytes','modifiedDatenum'});
end

function touch(path,mode)
fid=fopen(path,mode); assert(fid>=0); fwrite(fid,uint8(1)); fclose(fid);
end

function removeTree(path)
if isfolder(path), rmdir(path,'s'); end
end
