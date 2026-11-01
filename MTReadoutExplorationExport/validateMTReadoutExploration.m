function validateMTReadoutExploration()
%VALIDATEMTREADOUTEXPLORATION Synthetic four-task integration checks.
% Uses the real project selector; writes only within a temporary directory.
assert(exist('selectAnalysisFiles','file')==2,'Put the current project selector on the path.');
hadPreference=ispref('MTReadoutExploration','projectRoot'); oldRoot='';
if hadPreference, oldRoot=getpref('MTReadoutExploration','projectRoot'); end
restore=onCleanup(@()restorePreference(hadPreference,oldRoot)); %#ok<NASGU>
root=tempname; mkdir(root); cleanup=onCleanup(@()rmdir(root,'s')); %#ok<NASGU>
makeFixture(root,'IDS',10,true,true);
makeFixture(root,'IDQ',120,true,false);
makeFixture(root,'IDQ',120,false,false);
makeFixture(root,'IDJ',120,true,false);
makeFixture(root,'IDR',179,false,false);
makeFixture(root,'IDR',180,false,false);
out=fullfile(root,'export');
bundle=exportMTReadoutExploration('ProjectRoot',root,'OutputFolder',out,'IncludeAnalyses',false);
assert(isfile(bundle)); M=readtable(fullfile(out,'SourceManifest.csv'),'TextType','string');
assert(height(M)==5,'Expected IDS, randomized IDQ, IDJ, and two IDR probe files.');
assert(sum(M.n_trials)==9,'Top-up/fixed-IDQ exclusion or trial selection failed.');
for k=1:height(M)
  A=load(fullfile(out,char(M.bundle_file(k))),'source'); E=A.source;
  assert(ischar(E.meta.domain)&&~any(E.topup_trial));
  assert(isequal(size(E.noise_pc),[E.meta.n_trials 2 3 20]));
  assert(all(E.correct==(E.choice_side==E.change_side)));
  ip=E.bin_start_ms>=500 & E.bin_stop_ms<=750;
  w=reshape(E.bin_frame_count(ip)/sum(E.bin_frame_count(ip)),1,1,1,[]);
  actual=E.deterministic_pre_pc+sum(E.noise_pc(:,:,:,ip).*w,4);
  use=E.stream_mask;
  assert(max(abs(actual(use)-E.pre_window_total_pc(use)))<2e-5,'Prestep reconstruction differs.');
  if strcmp(E.meta.domain,'IDQ')
    assert(all(E.randomized_base_direction==1));
    assert(all(abs(double(E.directions_deg(:,1,1))-150)<1e-6));
    assert(all(abs(double(E.deterministic_step_pc(:,1,1))-22.071)<1e-4));
    assert(all(abs(double(E.noise_pc(:,1,2,1))-2)<1e-6),'IDQ positive probe order is wrong.');
  elseif strcmp(E.meta.domain,'IDJ')
    assert(E.meta.n_excluded_topup==1 && E.meta.n_trials==1);
    assert(abs(E.allocation_fraction(1,1)-.33330)<1e-8);
    assert(max(abs(E.nominal_allocation_fraction(1,:)-1/3))<1e-12);
  elseif strcmp(E.meta.domain,'IDR')
    if E.probe_offset_deg(1)==179
      assert(all(E.stream_mask(:,1,3))); assert(all(E.yoking_group(:,1,3)==2));
      assert(isequal(E.noise_pc(:,:,2,:),E.noise_pc(:,:,3,:)),'Yoking changed.');
    else
      assert(all(~E.stream_mask(:,1,3)),'180-degree IDR probe was duplicated.');
    end
  elseif strcmp(E.meta.domain,'IDS')
    assert(strcmp(E.meta.input_representation,'existing compact session analysis'));
    assert(isstruct(E.observed_columns),'MATLAB table was not serialized.');
  end
end
% A second export uses saved source shards and should not parse trials again.
exportMTReadoutExploration('ProjectRoot',root,'OutputFolder',out,'IncludeAnalyses',false);
M=readtable(fullfile(out,'SourceManifest.csv'),'TextType','string');
assert(all(M.cache_status=="reused"),'Unchanged source exports were not reused.');
fprintf('Four-domain export, IDQ ordering, IDR yoking, IDJ top-ups/thirds, and cache checks passed.\n');
end

function restorePreference(had,old)
if had
  setpref('MTReadoutExploration','projectRoot',old);
elseif ispref('MTReadoutExploration','projectRoot')
  rmpref('MTReadoutExploration','projectRoot');
end
end

function makeFixture(root,domain,offset,random,compact)
sid=sprintf('Fixture_%s_Meetz_%d_%d',domain,offset,random);
if any(strcmp(domain,{'IDS','IDR'}))
  folder=fullfile(root,domain,'Data',sprintf('Probe%d',offset),'ProbeSessions');
else
  folder=fullfile(root,domain,'Data','FullSessions');
end
if ~isfolder(folder), mkdir(folder); end
sessionHeader=struct('sourceDomain',domain,'animal','Meetz','fileName',sid, ...
  'frameRateHz',120,'preStepMS',750,'stepMS',250,'randomBaseDir',random,'prefDirDeg',30);
sessionProbeHeader=struct('sourceDomain',domain,'sourceSessionID',sid, ...
  'probeOffsetDeg',offset,'probeDirDeg',offset,'probeCohNoisePC',7.071,'randomBaseDir',random);
trials=cell(1,2); nf=120; N=zeros(2,3,nf,2);
for k=1:2
  D=struct('instructTrial',false,'randomizeDriftDir',true,'preStepMS',750,'stepMS',250, ...
    'stepSide',0,'changeSide',0,'baseCohPC',7.071,'stepCohPC',22.071, ...
    'cohNoisePC',7.071,'stepHasNoise',true,'stepSign',1,'probeOffsetDeg',offset, ...
    'driftDirDeg',150,'dirIndex',1,'dirDeg',150,'baseDirDeg',30,'randomBaseDir',random);
  if strcmp(domain,'IDQ'), D.stepCohPC=15; end
  if strcmp(domain,'IDR')
    D.baseCohPC=35; D.stepCohPC=55; D.changeIndex=1; D.cohNoise=true;
    D.probeDirDeg=offset; D.probeCohNoisePC=7.071;
  elseif strcmp(domain,'IDJ')
    D.topUpTrial=k==2; D.prestepCohPC=7.071; D.alloc0=.33330; D.alloc1=.33335; D.alloc2=.33335;
    D.totalStepCohPC=21; D.stepEndCoh0PC=7.071+21*D.alloc0;
    D.stepEndCoh1PC=7.071+21*D.alloc1; D.stepEndCoh2PC=7.071+21*D.alloc2;
    D.signalDir0Deg=45; D.signalDir1Deg=165; D.signalDir2Deg=285;
  end
  tr=struct('trial',struct('data',D),'trialCertify',struct('data',0), ...
    'extendedEOT',struct('data',0),'trialEnd',struct('data',0), ...
    'targetChosen',struct('data',0),'trialNumber',k);
  for prefix={'change','noChange'}
    tr.([prefix{1} 'TimesMS'])=struct('data',(0:50:950)');
    if strcmp(domain,'IDR'), roles={'Pref','Probe'};
    elseif strcmp(domain,'IDJ'), roles={'Noise0','Noise1','Noise2'};
    elseif strcmp(domain,'IDQ'), roles={'CohsPC0','CohsPC1','CohsPC2'};
    else, roles={'Drift','Probe0','Probe1'}; end
    for s=1:numel(roles)
      if strcmp(domain,'IDQ'), field=[prefix{1} roles{s}]; value=2*(s-2);
      else, field=[prefix{1} roles{s} 'CohsPC']; value=s-1; end
      tr.(field)=struct('data',repmat(value,20,1));
    end
  end
  trials{k}=tr;
  for s=1:3, N(:,s,:,k)=s-1; end
end
path=fullfile(folder,[sid '.mat']);
if any(strcmp(domain,{'IDS','IDR'}))
  probeTrials=trials; trialIdx=1:2; %#ok<NASGU>
  save(path,'sessionHeader','sessionProbeHeader','probeTrials','trialIdx','-v7');
else
  save(path,'sessionHeader','trials','-v7');
end
if compact
  analysisFolder=fullfile(fileparts(folder),'SessionAnalysis'); mkdir(analysisFolder);
  T=table((1:2)',ones(2,1),ones(2,1),true(2,1),repmat(7.071,2,1), ...
    repmat(22.071,2,1),true(2,1),repmat(150,2,1),repmat(offset,2,1), ...
    repmat(7.071,2,1),repmat(2,2,1),'VariableNames',{'parentTrialIdx', ...
    'chosenSideIndex','sideIndex','correct','baseCoh','stepEndpointCoh', ...
    'hasStepNoise','driftDirDeg','probeOffsetDeg','cohNoisePC','stepSignIndex'});
  sessionAnalysis=struct('sessionHeader',sessionHeader,'sessionProbeHeader',sessionProbeHeader, ...
    'trialTable',T,'noiseBySideStream',N); %#ok<NASGU>
  save(fullfile(analysisFolder,[sid '_sessionAnalysis.mat']), ...
    'sessionAnalysis','sessionHeader','sessionProbeHeader','-v7');
end
end
