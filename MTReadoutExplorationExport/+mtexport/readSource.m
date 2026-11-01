function E = readSource(item,binMS)
%READSOURCE Export stimulus inputs, never population-model predictions.
E=[]; compactNote='';
if ~isempty(item.compact)
  a=dir(item.compact); b=dir(item.source);
  if a.datenum>=b.datenum
    try
      S=load(item.compact,'sessionAnalysis');
      if isfield(S,'sessionAnalysis'), E=readCompact(S.sessionAnalysis,item); end
    catch ME
      compactNote=['Compact analysis not used: ' ME.message]; E=[];
    end
  end
end
if isempty(E), E=readRaw(item); E.meta.input_representation='source trial events';
else, E.meta.input_representation='existing compact session analysis'; end
if ~isempty(compactNote), E.meta.compact_fallback_note=compactNote; end
E=finish(E,binMS);
end

function E=emptyRecord(n,H,item,nFrames)
E=struct(); E.meta=struct('schema_version',1,'domain',item.domain, ...
  'animal',item.animal,'session_id',item.session,'side_index_base',1, ...
  'source_name',baseName(item.source),'coherence_units','percentage points', ...
  'direction_units','degrees','noise_convention','recorded zero-centered coherence fluctuation', ...
  'n_excluded_topup',0,'n_excluded_nonrandom_idq',0);
rate=double(mtexport.get(H,'frameRateHz',NaN));
pre=double(mtexport.get(H,'preStepMS',NaN)); step=double(mtexport.get(H,'stepMS',NaN));
assert(all(isfinite([rate pre step]))&&rate>0&&pre>=250&&step>0,'Missing/invalid timing metadata.');
E.frame_time_ms=(0:nFrames-1)*1000/rate;
E.meta.frame_rate_hz=rate; E.meta.step_start_ms=pre; E.meta.step_duration_ms=step;
E.meta.pre_window_ms=[pre-250 pre]; E.meta.step_window_ms=[pre pre+step];
E.parent_trial_index=zeros(n,1); E.source_trial_number=nan(n,1);
E.change_side=zeros(n,1,'uint8'); E.choice_side=zeros(n,1,'uint8');
E.correct=false(n,1); E.has_step_noise=false(n,1); E.is_increment=true(n,1);
E.randomized_base_direction=nan(n,1); E.probe_offset_deg=nan(n,1);
E.allocation_fraction=zeros(n,3); E.nominal_allocation_fraction=zeros(n,3);
E.directions_deg=nan(n,2,3,'single'); E.stream_mask=false(n,2,3);
E.yoking_group=zeros(n,2,3,'uint8');
E.deterministic_pre_pc=nan(n,2,3,'single'); E.deterministic_step_pc=nan(n,2,3,'single');
E.noise_frame_pc=nan(n,2,3,nFrames,'single'); E.topup_trial=false(n,1);
E.meta.topups_excluded=true;
end

function E=readCompact(A,item)
E=[];
if ~isfield(A,'trialTable')||~istable(A.trialTable), return; end
T=A.trialTable; H=A.sessionHeader; n=height(T); names=T.Properties.VariableNames;
if strcmp(item.domain,'IDQ')&&~isfinite(item.randomBase), return; end
if strcmp(item.domain,'IDJ')
  req={'sourceTrialIndex','chosenSideIndex','stepSideIndex','prestepCohPC', ...
    'stepEndCoh0PC','stepEndCoh1PC','stepEndCoh2PC','signalDir0Deg','signalDir1Deg', ...
    'signalDir2Deg','alloc0','alloc1','alloc2','topUpTrial'};
  current=isfield(A,'topUpTrialsExcluded')&&isequal(A.topUpTrialsExcluded,true) && ...
    all(ismember(req,names))&&all(~logical(T.topUpTrial));
  if ~current || ~all(isfield(A,{'fullChangeNoiseByFrameTrial','fullNoChangeNoiseByFrameTrial'})), return; end
  C=A.fullChangeNoiseByFrameTrial; U=A.fullNoChangeNoiseByFrameTrial; nf=size(C,1);
  E=emptyRecord(n,H,item,nf); E.meta.n_excluded_topup=double(mtexport.get(A,'nExcludedTopUpTrials',0));
  E.parent_trial_index=double(T.sourceTrialIndex); E.source_trial_number=double(T.sourceTrialNumber);
  E.change_side=uint8(T.stepSideIndex); E.choice_side=uint8(T.chosenSideIndex);
  E.correct=logical(T.correct); E.has_step_noise=logical(T.stepHasNoise);
  E.randomized_base_direction(:)=1; E.allocation_fraction=double(T{:,{'alloc0','alloc1','alloc2'}});
  for k=1:n
    E=put(E,k,double(T{k,{'signalDir0Deg','signalDir1Deg','signalDir2Deg'}}), ...
      repmat(double(T.prestepCohPC(k)),1,3), ...
      double(T{k,{'stepEndCoh0PC','stepEndCoh1PC','stepEndCoh2PC'}}), ...
      reshape(C(:,k,:),nf,3),reshape(U(:,k,:),nf,3),[1 2 3]);
  end
else
  req={'parentTrialIdx','chosenSideIndex','sideIndex','baseCoh','stepEndpointCoh', ...
    'hasStepNoise','driftDirDeg','probeOffsetDeg','cohNoisePC'};
  if ~all(ismember(req,names))||~isfield(A,'noiseBySideStream'), return; end
  N=A.noiseBySideStream; nf=size(N,3);
  assert(size(N,1)==2&&size(N,2)==3&&size(N,4)==n,'Unexpected compact noise-array dimensions.');
  E=emptyRecord(n,H,item,nf); E.parent_trial_index=double(T.parentTrialIdx);
  E.change_side=uint8(T.sideIndex); E.choice_side=uint8(T.chosenSideIndex);
  E.correct=logical(T.correct); E.has_step_noise=logical(T.hasStepNoise);
  if ismember('stepSignIndex',names), E.is_increment=double(T.stepSignIndex)==2; end
  E.probe_offset_deg=double(T.probeOffsetDeg); E.allocation_fraction(:,1)=1;
  E.randomized_base_direction(:)=1;
  if strcmp(item.domain,'IDQ')
    assert(isfinite(item.randomBase),'IDQ compact analysis lacks explicit randomized-base metadata; use source events.');
    E.randomized_base_direction(:)=item.randomBase;
  end
  for k=1:n
    dirs=mod(double(T.driftDirDeg(k))+[0 double(T.probeOffsetDeg(k)) -double(T.probeOffsetDeg(k))],360);
    base=[double(T.baseCoh(k)) repmat(double(T.cohNoisePC(k)),1,2)];
    post=base; post(1)=double(T.stepEndpointCoh(k));
    c=double(E.change_side(k));
    E=put(E,k,dirs,base,post,reshape(N(c,:,:,k),3,nf)', ...
      reshape(N(3-c,:,:,k),3,nf)',[1 2 3]);
  end
end
E.observed_columns=mtexport.plain(T);
if isfield(A,'configurationTable'), E.configuration=mtexport.plain(A.configurationTable); end
end

function E=readRaw(item)
vars=who('-file',item.source);
wanted=intersect(vars,{'sessionHeader','sessionProbeHeader','header','trials','probeTrials','trialIdx'},'stable');
S=load(item.source,wanted{:}); H=S.sessionHeader; P=struct();
if isfield(S,'sessionProbeHeader'), P=S.sessionProbeHeader; end
if isfield(S,'probeTrials'), trials=S.probeTrials; else, trials=S.trials; end
good=false(numel(trials),1); completed=good; nTop=0; nFixed=0;
for k=1:numel(trials)
  tr=trials{k};
  if ~isfield(tr,'trial')||~isfield(tr.trial,'data')||~isfield(tr,'trialCertify')||~isfield(tr,'targetChosen'), continue; end
  D=tr.trial.data; eot=outcome(tr,item.domain);
  if mtexport.get(tr,'trialCertify',1)~=0||~ismember(eot,[0 1]), continue; end
  completed(k)=true;
  if logical(mtexport.get(D,'instructTrial',false)), continue; end
  top=logical(mtexport.get(D,{'topUpTrial','isTopUpTrial','topUp','isTopUp'},false));
  if strcmp(item.domain,'IDJ'), assert(isfield(D,'topUpTrial'),'IDJ trial lacks required topUpTrial flag.'); end
  if top, nTop=nTop+1; continue; end
  if strcmp(item.domain,'IDQ')
    mode=mtexport.get(D,'randomBaseDir',mtexport.get(P,'randomBaseDir',mtexport.get(H,'randomBaseDir',false)));
    if ~logical(mode), nFixed=nFixed+1; continue; end
  end
  good(k)=true;
end
use=find(good); n=numel(use);
if n==0
  error('MTReadoutExport:NoEligibleTrials','No eligible completed non-top-up trials (or no randomized IDQ trials).');
end
rate=double(mtexport.get(H,'frameRateHz',NaN)); pre=double(mtexport.get(H,'preStepMS',NaN));
step=double(mtexport.get(H,'stepMS',NaN)); nf=round((pre+step)*rate/1000);
E=emptyRecord(n,H,item,nf); E.meta.n_excluded_topup=nTop; E.meta.n_excluded_nonrandom_idq=nFixed;
parent=(1:numel(trials))';
if isfield(S,'trialIdx')
  idx=double(S.trialIdx(:));
  if numel(idx)==numel(trials), parent=idx;
  elseif numel(idx)==nnz(completed)
    parent=nan(numel(trials),1); parent(completed)=idx;
  elseif numel(idx)==nnz(good)
    parent=nan(numel(trials),1); parent(good)=idx;
  else, error('Ambiguous trialIdx contract in %s.',item.source); end
end
E.parent_trial_index=parent(use);
for row=1:n
  tr=trials{use(row)}; D=tr.trial.data;
  E.source_trial_number(row)=double(mtexport.get(tr,'trialNumber',NaN));
  E.correct(row)=outcome(tr,item.domain)==0;
  E.choice_side(row)=uint8(double(mtexport.get(tr,'targetChosen',NaN))+1);
  if strcmp(item.domain,'IDJ')
    c=double(D.stepSide)+1; E.change_side(row)=uint8(c);
    E.has_step_noise(row)=logical(D.stepHasNoise); E.randomized_base_direction(row)=1;
    dirs=double([D.signalDir0Deg D.signalDir1Deg D.signalDir2Deg]);
    base=repmat(double(D.prestepCohPC),1,3);
    post=double([D.stepEndCoh0PC D.stepEndCoh1PC D.stepEndCoh2PC]);
    E.allocation_fraction(row,:)=double([D.alloc0 D.alloc1 D.alloc2]);
    fields={'Noise0','Noise1','Noise2'};
  elseif strcmp(item.domain,'IDR')
    E.change_side(row)=uint8(double(D.changeSide)+1);
    E.has_step_noise(row)=logical(mtexport.get(D,'cohNoise',true));
    E.is_increment(row)=double(mtexport.get(D,{'changeIndex','stepIndex'},1))==1;
    drift=mtexport.get(H,'prefDirDeg',NaN);
    if ~isfinite(drift)&&isfield(S,'header'), drift=mtexport.get(S.header,'prefDirDeg',NaN); end
    assert(isfinite(drift),'IDR lacks preferred-direction metadata.');
    offset=double(mtexport.get(D,'probeDirDeg',mtexport.get(P,'probeDirDeg',NaN)));
    assert(isfinite(offset),'IDR lacks probe offset.');
    pedestal=double(mtexport.get(D,'probeCohNoisePC',mtexport.get(P,'probeCohNoisePC',NaN)));
    E.randomized_base_direction(row)=0; E.allocation_fraction(row,1)=1;
    if offset<0 && ~E.has_step_noise(row)
      dirs=double(drift); base=double(D.baseCohPC); post=double(D.stepCohPC); fields={'Pref'};
    elseif abs(offset-180)<1e-6
      dirs=mod(double(drift)+[0 180],360); base=[double(D.baseCohPC) pedestal];
      post=[double(D.stepCohPC) pedestal]; fields={'Pref','Probe'};
    else
      dirs=mod(double(drift)+[0 offset -offset],360);
      base=[double(D.baseCohPC) pedestal pedestal]; post=[double(D.stepCohPC) pedestal pedestal];
      fields={'Pref','Probe','Probe'};
    end
    E.probe_offset_deg(row)=offset;
  else
    E.change_side(row)=uint8(double(mtexport.get(D,{'stepSide','changeSide'},NaN))+1);
    E.has_step_noise(row)=logical(D.stepHasNoise); E.allocation_fraction(row,1)=1;
    E.randomized_base_direction(row)=1;
    canonical=isfield(tr,'changeDriftCohsPC')&&isfield(D,'driftDirDeg')&&isfield(D,'baseCohPC');
    if canonical
      drift=double(D.driftDirDeg); offset=double(D.probeOffsetDeg);
      base=[double(D.baseCohPC) double(D.cohNoisePC) double(D.cohNoisePC)];
      post=base; post(1)=double(D.stepCohPC); fields={'Drift','Probe0','Probe1'};
      E.is_increment(row)=double(D.stepSign)==1;
    else
      assert(strcmp(item.domain,'IDQ'),'Unsupported IDS trial descriptor.');
      drift=double(D.dirDeg); offset=120; di=double(D.dirIndex);
      assert(ismember(di,0:2),'Invalid native IDQ dirIndex.');
      fields={sprintf('CohsPC%d',di),sprintf('CohsPC%d',mod(di+1,3)),sprintf('CohsPC%d',mod(di+2,3))};
      base=repmat(double(D.cohNoisePC),1,3); post=base; post(1)=base(1)+double(D.stepCohPC);
    end
    dirs=mod(drift+[0 offset -offset],360); E.probe_offset_deg(row)=offset;
  end
  assert(abs(double(mtexport.get(D,'preStepMS',pre))-pre)<1e-6 && ...
    abs(double(mtexport.get(D,'stepMS',step))-step)<1e-6,'Trial timing differs within source.');
  nativeIDQ=strcmp(item.domain,'IDQ')&&startsWith(fields{1},'CohsPC');
  change=rawTrace(tr,'change',fields,E.frame_time_ms,nativeIDQ,strcmp(item.domain,'IDJ'));
  unchanged=rawTrace(tr,'noChange',fields,E.frame_time_ms,nativeIDQ,strcmp(item.domain,'IDJ'));
  groups=1:numel(dirs);
  if strcmp(item.domain,'IDR')&&numel(dirs)==3, groups=[1 2 2]; end
  E=put(E,row,dirs,base,post,change,unchanged,groups);
end
end

function E=put(E,k,dirs,base,post,C,U,groups)
s=numel(dirs); nf=numel(E.frame_time_ms); c=double(E.change_side(k));
assert(ismember(c,[1 2])&&ismember(double(E.choice_side(k)),[1 2]),'Invalid side indices.');
assert(E.correct(k)==(E.choice_side(k)==E.change_side(k)),'Choice and outcome disagree.');
assert(isequal(size(C),[nf s])&&isequal(size(U),[nf s]),'Noise-array shape mismatch.');
for j=1:2
  E.directions_deg(k,j,1:s)=reshape(single(dirs),1,1,s); E.stream_mask(k,j,1:s)=true;
  E.yoking_group(k,j,1:s)=reshape(uint8(groups),1,1,s);
  E.deterministic_pre_pc(k,j,1:s)=reshape(single(base),1,1,s);
  E.deterministic_step_pc(k,j,1:s)=reshape(single(base),1,1,s);
end
E.deterministic_step_pc(k,c,1:s)=reshape(single(post),1,1,s);
E.noise_frame_pc(k,c,1:s,:)=reshape(single(C'),1,1,s,nf);
E.noise_frame_pc(k,3-c,1:s,:)=reshape(single(U'),1,1,s,nf);
end

function M=rawTrace(tr,prefix,roles,t,nativeIDQ,previousEvent)
times=double(mtexport.value(tr.([prefix 'TimesMS']))); times=times(:);
M=nan(numel(t),numel(roles),'single'); frameMS=t(2)-t(1);
for s=1:numel(roles)
  if nativeIDQ, field=[prefix roles{s}]; else, field=[prefix roles{s} 'CohsPC']; end
  v=double(mtexport.value(tr.(field))); v=v(:);
  assert(numel(v)==numel(times),'Noise/time event count mismatch.');
  if previousEvent
    for f=1:numel(t)
      q=find(times<=t(f)+1e-6,1,'last'); assert(~isempty(q),'Missing initial noise event.'); M(f,s)=single(v(q));
    end
  else
    for q=1:numel(times)
      a=max(1,floor(times(q)/frameMS)+1);
      if q<numel(times), b=floor(times(q+1)/frameMS)+1; else, b=numel(t)+1; end
      b=min(b,numel(t)+1);
      if a<=numel(t)&&b>a, M(a:b-1,s)=single(v(q)); end
    end
  end
end
assert(all(isfinite(M),'all'),'Noise events do not cover the full trial.');
end

function E=finish(E,binMS)
n=numel(E.parent_trial_index); assert(n>0,'Empty export source.');
assert(all(isfinite(E.parent_trial_index))&&numel(unique(E.parent_trial_index))==n,'Nonunique/invalid parent indices.');
if strcmp(E.meta.domain,'IDQ'), assert(all(E.randomized_base_direction==1),'Fixed/unknown IDQ entered export.'); end
E.meta.n_trials=n; E.meta.n_fit_trials=nnz(E.has_step_noise&E.is_increment);
E.nominal_allocation_fraction=E.allocation_fraction;
if strcmp(E.meta.domain,'IDJ')
  assert(all(abs(sum(E.allocation_fraction,2)-1)<1e-4),'IDJ allocations do not sum to 1.');
  equal=max(abs(E.allocation_fraction-1/3),[],2)<1e-4;
  E.nominal_allocation_fraction(equal,:)=1/3;
end
pre=E.meta.step_start_ms; stop=pre+E.meta.step_duration_ms; t=E.frame_time_ms;
ip=t>=pre-250 & t<pre; is=t>=pre & t<stop;
assert(nnz(ip)==nnz(is),'Unequal pre/post frame counts.');
E.pre_window_total_pc=E.deterministic_pre_pc+mean(E.noise_frame_pc(:,:,:,ip),4);
E.step_window_total_pc=E.deterministic_step_pc+mean(E.noise_frame_pc(:,:,:,is),4);
if binMS==0
  edges=[t stop];
else
  edges=unique([0:binMS:stop pre-250 pre stop]);
end
nb=numel(edges)-1; E.noise_pc=nan(n,2,3,nb,'single'); counts=zeros(1,nb);
for k=1:nb
  use=t>=edges(k) & t<edges(k+1); counts(k)=nnz(use);
  assert(any(use),'Requested bin width is below the frame resolution.');
  E.noise_pc(:,:,:,k)=mean(E.noise_frame_pc(:,:,:,use),4);
end
E.bin_start_ms=edges(1:end-1); E.bin_stop_ms=edges(2:end);
E.bin_frame_count=counts; E.meta.bin_ms=binMS;
E.meta.temporal_representation='mean recorded coherence noise per bin; deterministic pre/post levels stored separately';
E.meta.yoking_contract='Group IDs apply within each RDK. They do not imply cross-RDK correlation.';
E.meta.stream_convention='drift,+probe,-probe for IDR/IDQ/IDS; recorded signal0,signal1,signal2 for IDJ; inactive streams masked';
E.meta.has_negative_window_coherence=any(E.pre_window_total_pc(:)<-1e-4)||any(E.step_window_total_pc(:)<-1e-4);
E.meta.min_window_coherence_pc=min([E.pre_window_total_pc(E.stream_mask);E.step_window_total_pc(E.stream_mask)]);
E.meta.max_window_coherence_budget_pc=max([sum(E.pre_window_total_pc,3,'omitnan'); ...
  sum(E.step_window_total_pc,3,'omitnan')],[],'all');
E.meta.coherence_audit_note='Recorded signed values are preserved without clipping or dot-budget repair.';
E=rmfield(E,{'noise_frame_pc','frame_time_ms'});
end

function eot=outcome(tr,domain)
if strcmp(domain,'IDR'), names={'trialEnd','extendedEOT'};
else, names={'extendedEOT','trialEnd'}; end
eot=mtexport.get(tr,names,NaN);
end

function name=baseName(path)
[~,stem,ext]=fileparts(path); name=[stem ext];
end
