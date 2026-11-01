function validateMTReadoutPNormLocal()
%VALIDATEMTREADOUTPNORMLOCAL Small integration checks using the real selector.
% No real recordings are modified; fixtures and outputs live in tempname.
validateMTReadoutPNorm;
assert(exist('selectAnalysisFiles','file')==2,'Add the project selector to the path.');
root=tempname; mkdir(root); cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
for offset=[0 120 179 180]
  fixture(root,offset,sprintf('FixtureMeetz_Probe%d.mat',offset),'Meetz','IDS');
end
fixture(root,120,'FixtureNeesha_Probe120.mat','Neesha','IDS');
fixture(root,120,'FixtureIDQ_Probe120.mat','Meetz','IDQ');
fixture(root,120,'ProbeSessionIDReadout_Meetz_20260114_Probe120.mat','Meetz','IDS');
[files,info]=mtpnorm.selectLocalFiles(root,"Meetz",false);
assert(numel(files)==4,'Selector filters/exclusions were not preserved.');
assert(any(info.probeDirDeg==179) && any(info.probeDirDeg==180), ...
  '179/180 offsets were unexpectedly combined or omitted.');
cfg=mtpnorm.defaultConfig(); cfg.trialBatchSize=1;
cacheFolder=fullfile(root,'cache');
[C,M]=mtpnorm.prepareLocalCache(files,cacheFolder,cfg,false);
assert(all(M.cacheStatus=="extracted"));
[C2,M]=mtpnorm.prepareLocalCache(files,cacheFolder,cfg,false);
assert(all(M.cacheStatus=="reused"));
assert(isequaln(C.preCohPC,C2.preCohPC));
marker=pi; save(files{1},'marker','-append');
[C,M]=mtpnorm.prepareLocalCache(files,cacheFolder,cfg,false);
assert(nnz(M.cacheStatus=="extracted")==1,'Changed source was not re-extracted once.');
assert(~isfield(C,'rawStepNoisePC'),'Compact cache retained unused frame traces.');
R=mtpnorm.fullResponses(C,cfg); Z=mtpnorm.noiseFields(cfg,16);
[pr,a]=mtpnorm.predict(R,5,12,Z);
[batched,b]=mtpnorm.predictCache(C,cfg,5,12,Z);
assert(max(abs(pr-batched))<1e-12,'Batched prediction differs.');
assert(abs(a.negativeRateFraction-b.negativeRateFraction)<1e-12,'Clipping audit differs.');
idx=[2 4]; sub=mtpnorm.predictCache(C,cfg,5,12,Z,idx);
assert(max(abs(sub-pr(idx)))<1e-12,'Batch subset was ignored.');
result=struct('config',cfg,'trials',C.trials,'fitMask',C.trials.hasStepNoise & C.trials.isIncrement, ...
  'fits',table(5,12,'VariableNames',{'p','noiseSD'}), ...
  'gainStatus',"test only",'gainNotes',strings(0,1));
pc=pr; q=C.trials.changeRDK==2; pc(q)=1-pr(q);
result.predictions={struct('pRDK1',pr,'pCorrect',pc)};
out=fullfile(root,'out'); mkdir(out);
share=exportMTReadoutPNormSummary(result,out); assert(isfile(share));
A=load(fullfile(out,'Share','MTReadoutPNormSummary.mat'),'summary');
assert(~any(isfield(A.summary,{'trials','predictions','preCohPC','stepCohPC'})), ...
  'Trial-level data leaked into compact summary.');
assert(height(A.summary.coverage)==4);
fprintf('Local selector, incremental cache, batching, subset, and share-export checks passed.\n');
end

function fixture(root,offset,name,animal,domain)
folder=fullfile(root,'Data',sprintf('Probe%d',offset),'ProbeSessions');
if ~isfolder(folder), mkdir(folder); end
sid=sprintf('Fixture_%s_%s_%d',animal,domain,offset);
sessionHeader=struct('animal',animal,'sourceDomain',domain,'fileName',sid, ...
  'createdBy','convertIDSData','frameRateHz',120,'preStepMS',750,'stepMS',250, ...
  'nProbeDirections',1,'numberOfTrials',2);
sessionProbeHeader=struct('sourceDomain',domain,'sourceSessionID',sid, ...
  'probeOffsetDeg',offset,'probeTag',sprintf('Probe%d',offset), ...
  'parentNProbeDirections',1,'nNoiseTrials',2);
probeTrials=cell(1,2); trialIdx=[1 2];
for k=1:2
  td=struct('instructTrial',false,'randomizeDriftDir',true,'preStepMS',750,'stepMS',250, ...
    'stepSide',k-1,'baseCohPC',7.071,'stepCohPC',22.071,'cohNoisePC',7.071, ...
    'probeOffsetDeg',offset,'driftDirDeg',90*k,'stepHasNoise',true,'stepSign',1);
  tr=struct('trial',struct('data',td),'trialCertify',struct('data',0), ...
    'extendedEOT',struct('data',0),'targetChosen',struct('data',k-1));
  for prefix={'change','noChange'}
    tr.([prefix{1} 'TimesMS'])=struct('data',(0:50:950)');
    for role={'Drift','Probe0','Probe1'}
      tr.([prefix{1} role{1} 'CohsPC'])=struct('data',zeros(20,1));
    end
  end
  probeTrials{k}=tr;
end
save(fullfile(folder,name),'sessionHeader','sessionProbeHeader','probeTrials','trialIdx','-v7');
end
