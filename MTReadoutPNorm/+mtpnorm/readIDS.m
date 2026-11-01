function cache = readIDS(sourceFiles,cfg)
%READIDS Extract corrected IDS raw sessions or ProbeSessions.
% Supports top-level trials (raw converted IDS) or probeTrials (ProbeSession).
% Recorded coherence traces are zero-centered imposed noise, not total
% coherence. Probe pedestal is cohNoisePC, added exactly once. Corrected
% baseCohPC and stepCohPC in trial.data are authoritative.
sourceFiles=cellstr(string(sourceFiles(:)));
rows=cell(0,14); pre=cell(0,2); post=cell(0,2); dirs=cell(0,1);
rawStep=cell(0,2); auditRows=cell(0,8); seen=strings(0,1);
for f=1:numel(sourceFiles)
  path=sourceFiles{f};
  % Avoid loading unrelated top-level data from the large source files.
  names=who('-file',path);
  wanted=intersect(names,{'sessionHeader','sessionProbeHeader','trials','probeTrials','trialIdx'},'stable');
  S=load(path,wanted{:}); H=S.sessionHeader;
  if isfield(H,'animal'), animal=string(H.animal); else, animal=cfg.animal; end
  if ~strcmpi(animal,cfg.animal), continue; end
  if isfield(S,'probeTrials')
    ts=S.probeTrials; P=S.sessionProbeHeader;
    if isfield(P,'sourceDomain')
      assert(strcmpi(string(P.sourceDomain),'IDS'),'Expected IDS, not IDQ/IDR.');
    end
    sid=string(fieldOr(P,'sourceSessionID',fieldOr(P,'sessionID',"")));
    parent=double(S.trialIdx(:));
  else
    assert(isfield(S,'trials'),'No trials or probeTrials in %s.',path);
    if isfield(H,'createdBy')
      assert(strcmpi(string(H.createdBy),'convertIDSData'), ...
        'Expected a convertIDSData recording, got %s.',string(H.createdBy));
    end
    ts=S.trials; sid=string(fieldOr(H,'fileName',"")); parent=(1:numel(ts))';
  end
  assert(iscell(ts),'Expected cell-array trials.');
  if strlength(sid)==0, [~,name]=fileparts(path); sid=string(name); end
  valid=false(numel(ts),1);
  for k=1:numel(ts), valid(k)=validTrial(ts{k}); end
  use=find(valid);
  if numel(parent)==numel(ts), parent=parent(use);
  else, assert(numel(parent)==numel(use),'Ambiguous ProbeSession trialIdx mapping.'); end
  nAdded=0; sourceOffsets=[]; correction=isfield(H,'dataCorrections');
  for z=1:numel(use)
    k=use(z); tr=ts{k}; td=tr.trial.data;
    key=sid+"/"+string(parent(z));
    assert(~ismember(key,seen),'Duplicate source trial %s: overlapping input files.',key);
    seen(end+1,1)=key; %#ok<AGROW>
    assert(logical(td.randomizeDriftDir),'Initial IDS model requires randomized drift.');
    preMS=double(td.preStepMS); stepMS=double(td.stepMS);
    assert(all(abs(cfg.preWindowMS-[preMS-250 preMS])<1e-6) && ...
      all(abs(cfg.stepWindowMS-[preMS preMS+stepMS])<1e-6) && abs(stepMS-250)<1e-6, ...
      'Timing differs from the agreed 500-750/750-1000-ms windows.');
    frameRate=double(value(H.frameRateHz)); frameMS=1000/frameRate;
    nFrames=round((preMS+stepMS)/frameMS); t=(0:nFrames-1)'*frameMS;
    ip=t>=cfg.preWindowMS(1) & t<cfg.preWindowMS(2);
    is=t>=cfg.stepWindowMS(1) & t<cfg.stepWindowMS(2);
    assert(nnz(ip)==nnz(is),'Pre/post windows have unequal frame counts.');
    c=double(td.stepSide)+1; u=3-c; choice=double(value(tr.targetChosen))+1;
    assert(ismember(c,[1 2])&&ismember(choice,[1 2]),'Invalid side coding.');
    correct=double(value(tr.extendedEOT))==0;
    assert(correct==(choice==c),'Choice and correctness disagree in %s.',key);
    base=double(td.baseCohPC); endpoint=double(td.stepCohPC);
    pedestal=double(td.cohNoisePC); offset=double(td.probeOffsetDeg);
    drift=mod(double(td.driftDirDeg),360); out=size(rows,1)+1;
    change=traceMatrix(tr,'change',nFrames,frameMS);
    unchanged=traceMatrix(tr,'noChange',nFrames,frameMS);
    bySide=cell(1,2); bySide{c}=change; bySide{u}=unchanged;
    for j=1:2
      noise=bySide{j}; before=mean(noise(ip,:),1)'; after=mean(noise(is,:),1)';
      before=before+[base;pedestal;pedestal];
      deterministic=[base;pedestal;pedestal];
      if j==c, deterministic(1)=endpoint; end
      after=after+deterministic;
      pre{out,j}=before; post{out,j}=after; %#ok<AGROW>
      if isfield(cfg,'keepRawStepNoise') && cfg.keepRawStepNoise
        rawStep{out,j}=noise(is,:); %#ok<AGROW>
      end
    end
    dirs{out,1}=mod([drift drift+offset drift-offset],360); %#ok<AGROW>
    rows(out,:)={sid,parent(z),string(path),animal,offset,drift,base,endpoint, ...
      abs(endpoint-base),c,choice,correct,logical(td.stepHasNoise),double(td.stepSign)==1}; %#ok<AGROW>
    nAdded=nAdded+1; sourceOffsets(end+1)=offset; %#ok<AGROW>
  end
  auditRows(end+1,:)={string(path),sid,numel(ts),nAdded, ...
    string(mat2str(unique(sourceOffsets))),correction,cfg.preWindowMS(1),cfg.preWindowMS(2)}; %#ok<AGROW>
end
assert(~isempty(rows),'No valid %s IDS trials.',cfg.animal);
T=cell2table(rows,'VariableNames',{'sessionID','parentTrialIdx','sourceFile','animal', ...
  'probeOffsetDeg','driftDirDeg','baseCohPC','stepEndpointCohPC','stepMagnitudePC', ...
  'changeRDK','choiceRDK','correct','hasStepNoise','isIncrement'});
for name={'correct','hasStepNoise','isIncrement'}, T.(name{1})=logical(T.(name{1})); end
sourceAudit=cell2table(auditRows,'VariableNames',{'sourceFile','sessionID', ...
  'candidateTrials','validTrials','offsetsDeg','hasCorrectionMetadata','preStartMS','preEndMS'});
cache=struct('schemaVersion',"MTReadoutPNormTrialCache-v2",'trials',T, ...
  'preCohPC',{pre},'stepCohPC',{post},'streamDirDeg',{dirs}, ...
  'rawStepNoisePC',{rawStep},'sourceAudit',sourceAudit,'config',cfg);
end

function good=validTrial(tr)
good=isstruct(tr)&&isfield(tr,'trial')&&isfield(tr.trial,'data')&& ...
  isfield(tr,'trialCertify')&&isfield(tr,'extendedEOT')&&isfield(tr,'targetChosen');
if ~good, return; end
td=tr.trial.data;
good=double(value(tr.trialCertify))==0 && ismember(double(value(tr.extendedEOT)),[0 1]) && ...
  ~logical(fieldOr(td,'instructTrial',false));
for name={'topUpTrial','isTopUpTrial','topUp','isTopUp','topupTrial','isTopupTrial'}
  if isfield(td,name{1}) && logical(value(td.(name{1}))), good=false; end
end
end

function M=traceMatrix(tr,prefix,nFrames,frameMS)
times=double(value(tr.([prefix 'TimesMS']))); times=times(:);
roles={'Drift','Probe0','Probe1'}; M=nan(nFrames,3);
for stream=1:3
  cohs=double(value(tr.([prefix roles{stream} 'CohsPC']))); cohs=cohs(:);
  assert(numel(times)==numel(cohs),'Trace/time length mismatch.');
  for q=1:numel(times)
    first=max(1,floor(times(q)/frameMS)+1);
    if q<numel(times), next=floor(times(q+1)/frameMS)+1; else, next=nFrames+1; end
    next=min(next,nFrames+1);
    if first<=nFrames && next>first, M(first:next-1,stream)=cohs(q); end
  end
end
assert(all(isfinite(M),'all'),'Coherence traces do not cover both analysis windows.');
end

function x=value(x)
while isstruct(x)&&isfield(x,'data'), x=x.data; end
end
function x=fieldOr(S,name,fallback)
if isfield(S,name), x=S.(name); else, x=fallback; end
end
