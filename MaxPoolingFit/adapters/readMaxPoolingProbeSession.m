function fragment = readMaxPoolingProbeSession(sourcePath)
%READMAXPOOLINGPROBESESSION Read one IDR, IDQ, or IDS ProbeSession.
%
% This adapter uses the established project conventions:
%   - stored *CohsPC traces are zero-centered imposed coherence noise;
%   - drift base/endpoint coherence comes from the trial descriptor;
%   - each probe has a positive pedestal equal to the trial/session noise
%     amplitude, added exactly once;
%   - IDR's one probe trace is duplicated into separate +offset and -offset
%     streams because the two displayed streams are yoked, not summed;
%   - normalized IDQ and corrected IDS trial descriptors are authoritative.

extractorVersion="ProbeSessionAdapter-v4-full-prestep";
S=load(sourcePath,'sessionHeader','sessionProbeHeader','probeTrials','trialIdx');
required={'sessionHeader','sessionProbeHeader','probeTrials','trialIdx'};
assert(all(isfield(S,required)),'ProbeSession lacks required variables: %s',sourcePath);
H=S.sessionHeader; P=S.sessionProbeHeader; trials=S.probeTrials;
domain=identifyDomain(H,P);
animal=identifyAnimal(H,sourcePath);
sessionID=string(fieldOr(P,'sourceSessionID',fieldOr(P,'sessionID',erase(string(sourcePath),".mat"))));

frameRate=scalarField(H,'frameRateHz');
preStepMS=scalarField(H,'preStepMS'); stepMS=scalarField(H,'stepMS');
frameMS=1000/frameRate; nFrames=round((preStepMS+stepMS)/frameMS);
tMS=(0:nFrames-1)'*frameMS;
stepFrames=find(tMS>=preStepMS & tMS<preStepMS+stepMS);
preFrames=find(tMS>=0 & tMS<preStepMS);

rows=cell(0,18); dirs=cell(0,2); integrated=cell(0,2);
preIntegrated=cell(0,2); drift=zeros(0,2); pre=zeros(0,2);
rawStep=cell(0,2); driver=cell(0,1);
rawNoiseValues=[]; validLocalIdx=zeros(0,1); out=0;
for k=1:numel(trials)
  tr=trials{k};
  [use,correct]=validCompletedTrial(tr,domain);
  if ~use, continue; end
  td=tr.trial.data; out=out+1; validLocalIdx(out,1)=k; %#ok<AGROW>
  [stepSignIndex,changeRDK,hasNoise,basePC,endpointPC,driftDeg,offsetDeg,pedestalPC]= ...
    trialMetadata(td,H,P,domain);
  choiceRDK=uint8(scalarEvent(tr,'targetChosen')+1);
  [changeNoise,noChangeNoise,streamDirs]=extractNoiseTraces( ...
    tr,domain,driftDeg,offsetDeg,nFrames,frameMS);
  changeStepNoise=changeNoise(stepFrames,:);
  noChangeStepNoise=noChangeNoise(stepFrames,:);
  rawNoiseValues=[rawNoiseValues; changeStepNoise(:); noChangeStepNoise(:)]; %#ok<AGROW>

  changeTotal=addDeterministicCoherence(changeNoise,endpointPC,pedestalPC);
  noChangeTotal=addDeterministicCoherence(noChangeNoise,basePC,pedestalPC);
  changePreTotal=addDeterministicCoherence(changeNoise,basePC,pedestalPC);
  noChangePreTotal=addDeterministicCoherence(noChangeNoise,basePC,pedestalPC);
  if changeRDK==1
    bySide={changeTotal(stepFrames,:),noChangeTotal(stepFrames,:)};
    preBySide={changePreTotal(preFrames,:),noChangePreTotal(preFrames,:)};
  else
    bySide={noChangeTotal(stepFrames,:),changeTotal(stepFrames,:)};
    preBySide={noChangePreTotal(preFrames,:),changePreTotal(preFrames,:)};
  end
  for j=1:2
    dirs{out,j}=streamDirs; %#ok<AGROW>
    rawStep{out,j}=bySide{j}; %#ok<AGROW>
    integrated{out,j}=mean(bySide{j},1,'omitnan')'; %#ok<AGROW>
    preIntegrated{out,j}=mean(preBySide{j},1,'omitnan')'; %#ok<AGROW>
  end
  % Temporal driver is zero-centered and aligned so positive favors the
  % changed RDK. Decrement trials are sign-reversed but retained only for
  % possible later analyses; primary kernel estimation selects increments.
  signStep=2*(stepSignIndex==2)-1;
  driver{out,1}=signStep.*sum(changeNoise(stepFrames,:)-noChangeNoise(stepFrames,:),2,'omitnan'); %#ok<AGROW>
  drift(out,:)=[driftDeg driftDeg]; pre(out,:)=[basePC basePC]; %#ok<AGROW>
  rows(out,:)={animal,domain,sessionID,NaN,choiceRDK,uint8(changeRDK), ...
    logical(correct),true,logical(hasNoise),uint8(stepSignIndex),stepSignIndex==2, ...
    animal+"_"+kernelDomainFamily(domain),basePC,endpointPC,abs(endpointPC-basePC), ...
    offsetDeg,driftDeg,string(sourcePath)}; %#ok<AGROW>
end

[parentTrialIdx,trialIdxContract]=resolveTrialIdx( ...
  S.trialIdx,validLocalIdx,numel(trials),sourcePath);
rows(:,4)=num2cell(parentTrialIdx);

T=cell2table(rows,'VariableNames',{'animal','domain','sessionID','trialID', ...
  'choiceRDK','changeRDK','correct','valid','hasStepNoise','stepSignIndex', ...
  'isIncrement','kernelFamily','baseCohPC','stepEndpointCohPC', ...
  'stepMagnitudePC','probeOffsetDeg','driftDirectionDeg','sourceFile'});
for f={'animal','domain','sessionID','kernelFamily','sourceFile'}
  T.(f{1})=string(T.(f{1}));
end
for f={'correct','valid','hasStepNoise','isIncrement'}
  T.(f{1})=logical(T.(f{1}));
end
fragment=struct('trials',T,'streamDirDeg',{dirs},'stepCohPC',{integrated}, ...
  'preCohPC',{preIntegrated},'driftDirDeg',drift,'preDriftCohPC',pre);
fragment.extractorVersion=extractorVersion;
fragment.raw=struct('stepStreamCohPC',{rawStep},'kernelDriverByFrame',{driver});
fragment.extractionAudit=struct('sourcePath',string(sourcePath),'domain',domain, ...
  'extractorVersion',extractorVersion, ...
  'nCandidateTrials',numel(trials),'nValidTrials',out, ...
  'rawNoiseMinPC',min(rawNoiseValues,[],'omitnan'), ...
  'rawNoiseMaxPC',max(rawNoiseValues,[],'omitnan'), ...
  'rawNoiseMeanPC',mean(rawNoiseValues,'omitnan'), ...
  'trialIdxContract',trialIdxContract, ...
  'traceConvention','stored trace treated as zero-centered imposed noise');
end

function [parentIdx,contract]=resolveTrialIdx(trialIdx,validLocalIdx,nCandidate,sourcePath)
% Support both historical ProbeSession trialIdx conventions.
trialIdx=double(trialIdx(:)); nValid=numel(validLocalIdx);
if numel(trialIdx)==nCandidate
  parentIdx=trialIdx(validLocalIdx);
  contract="one trialIdx per probeTrials entry";
elseif numel(trialIdx)==nValid
  parentIdx=trialIdx;
  contract="one trialIdx per valid analysis trial";
else
  error('readMaxPoolingProbeSession:TrialIdxLengthMismatch', ...
    ['Cannot align trialIdx in %s: %d trialIdx values, %d probeTrials ' ...
     'entries, and %d valid analysis trials.'], ...
    sourcePath,numel(trialIdx),nCandidate,nValid);
end
assert(all(isfinite(parentIdx))&&numel(unique(parentIdx))==numel(parentIdx), ...
  'Resolved parent trial indices must be finite and unique.');
end

function domain=identifyDomain(H,P)
if isfield(P,'sourceDomain'), domain=upper(string(P.sourceDomain));
elseif isfield(H,'sourceDomain'), domain=upper(string(H.sourceDomain));
elseif isfield(P,'probeOffsetDeg'), domain="IDS";
else, domain="IDR";
end
assert(ismember(domain,["IDR" "IDQ" "IDS"]),'Unrecognized source domain.');
end

function animal=identifyAnimal(H,path)
if isfield(H,'animal')&&strlength(string(H.animal))>0&&string(H.animal)~="unknown"
  animal=string(H.animal);
else
  names=["Meetz" "Neesha"];
  hit=arrayfun(@(name)contains(string(path),name),names);
  assert(sum(hit)==1,'Cannot identify animal from %s.',path); animal=names(hit);
end
end

function [use,correct]=validCompletedTrial(tr,domain)
use=isfield(tr,'trial')&&isfield(tr.trial,'data')&&isfield(tr,'trialCertify')&& ...
  scalarEvent(tr,'trialCertify')==0;
if ~use, correct=false; return; end
if domain=="IDR"
  if isfield(tr,'trialEnd'), eot=scalarEvent(tr,'trialEnd');
  elseif isfield(tr,'extendedEOT'), eot=scalarEvent(tr,'extendedEOT');
  else, use=false; correct=false; return; end
  instruction=false;
else
  if ~isfield(tr,'extendedEOT'), use=false; correct=false; return; end
  eot=scalarEvent(tr,'extendedEOT');
  instruction=isfield(tr.trial.data,'instructTrial')&&logical(tr.trial.data.instructTrial);
end
use=use&&~instruction&&ismember(eot,[0 1]); correct=eot==0;
end

function [signIndex,changeRDK,hasNoise,basePC,endpointPC,driftDeg,offsetDeg,pedestalPC]= ...
    trialMetadata(td,H,P,domain)
if domain=="IDR"
  if isfield(td,'changeIndex'), signIndex=double(td.changeIndex)+1;
  else, signIndex=double(td.stepIndex)+1; end
  changeRDK=double(td.changeSide)+1;
  hasNoise=~isfield(td,'cohNoise')||logical(td.cohNoise);
  basePC=double(td.baseCohPC); endpointPC=double(td.stepCohPC);
  driftDeg=scalarField(H,'prefDirDeg'); offsetDeg=double(P.probeDirDeg);
  pedestalPC=double(P.probeCohNoisePC);
else
  signIndex=double(td.stepSign)+1; changeRDK=double(td.stepSide)+1;
  hasNoise=logical(td.stepHasNoise); basePC=double(td.baseCohPC);
  endpointPC=double(td.stepCohPC); driftDeg=mod(double(td.driftDirDeg),360);
  offsetDeg=double(td.probeOffsetDeg); pedestalPC=double(td.cohNoisePC);
end
assert(ismember(signIndex,[1 2])&&ismember(changeRDK,[1 2]),'Invalid trial indices.');
end

function [changeNoise,noChangeNoise,dirs]=extractNoiseTraces(tr,domain,driftDeg,offsetDeg,nFrames,frameMS)
if domain=="IDR"
  changeDrift=fillTrace(tr,'changePrefCohsPC','changeTimesMS',nFrames,frameMS);
  noChangeDrift=fillTrace(tr,'noChangePrefCohsPC','noChangeTimesMS',nFrames,frameMS);
  changeProbe=fillTrace(tr,'changeProbeCohsPC','changeTimesMS',nFrames,frameMS);
  noChangeProbe=fillTrace(tr,'noChangeProbeCohsPC','noChangeTimesMS',nFrames,frameMS);
  if abs(offsetDeg-180)<1e-9
    changeNoise=[changeDrift changeProbe]; noChangeNoise=[noChangeDrift noChangeProbe];
    dirs=mod([driftDeg driftDeg+180],360);
  else
    changeNoise=[changeDrift changeProbe changeProbe];
    noChangeNoise=[noChangeDrift noChangeProbe noChangeProbe];
    dirs=mod([driftDeg driftDeg+offsetDeg driftDeg-offsetDeg],360);
  end
else
  roles={'Drift','Probe0','Probe1'};
  changeNoise=nan(nFrames,3); noChangeNoise=nan(nFrames,3);
  for r=1:3
    changeNoise(:,r)=fillTrace(tr,['change' roles{r} 'CohsPC'],'changeTimesMS',nFrames,frameMS);
    noChangeNoise(:,r)=fillTrace(tr,['noChange' roles{r} 'CohsPC'],'noChangeTimesMS',nFrames,frameMS);
  end
  dirs=mod([driftDeg driftDeg+offsetDeg driftDeg-offsetDeg],360);
end
end

function total=addDeterministicCoherence(noise,driftCohPC,pedestalPC)
total=noise+pedestalPC;
total(:,1)=noise(:,1)+driftCohPC;
end

function v=fillTrace(tr,cohField,timeField,nFrames,frameMS)
assert(isfield(tr,cohField)&&isfield(tr,timeField),'Missing %s or %s.',cohField,timeField);
cohs=dataValue(tr.(cohField)); cohs=double(cohs(:));
times=dataValue(tr.(timeField)); times=double(times(:));
n=min(numel(cohs),numel(times)); cohs=cohs(1:n); times=times(1:n); v=nan(nFrames,1);
for q=1:n
  first=max(1,floor(times(q)/frameMS)+1);
  if q<n, next=floor(times(q+1)/frameMS)+1; else, next=nFrames+1; end
  next=min(next,nFrames+1);
  if first<=nFrames&&next>first, v(first:next-1)=cohs(q); end
end
assert(all(isfinite(v)),'Trace %s does not cover all analysis frames.',cohField);
end

function value=scalarField(S,name)
assert(isfield(S,name),'Missing metadata field %s.',name); value=dataValue(S.(name));
if isnumeric(value)&&~isscalar(value), value=value(1); end
value=double(value);
end

function value=scalarEvent(tr,name)
value=dataValue(tr.(name)); value=double(value(1));
end

function value=dataValue(value)
while isstruct(value)&&isfield(value,'data'), value=value.data; end
end

function value=fieldOr(S,name,fallback)
if isfield(S,name), value=S.(name); else, value=fallback; end
end

function family=kernelDomainFamily(domain)
if domain=="IDR", family="IDR"; else, family="IDQIDS"; end
end
