function result = pilotIDLPNormProfile(varargin)
% pilotIDLPNormProfile  Pooled IDL p-norm estimation and uncertainty.
% TwoAllocationMode=false selects the original style; true selects 100:0/50:50.
%
% R = pilotIDLPNormProfile('Animal','Meetz');
% R = pilotIDLPNormProfile('Animal','Neesha','NBootstrap',1000);
% R = pilotIDLPNormProfile('Animal','Meetz','IncludeSession100',true);
% R = pilotIDLPNormProfile('Animal','Neesha','LapseMode','session', ...
%     'NBootstrap',100);
% R = pilotIDLPNormProfile('Animal','Neesha','LapseMode','fixed', ...
%     'LapseRate',.10,'NBootstrap',0);
%
% q = 1/p: q=0 is max, q=.5 is quadrature, q=1 is linear.
% PRIMARY (JointAll): for every q, fit one sigma per session using ALL
% allocations, then sum those session LLs. The pooled peak estimates one q
% per animal. Pedestals, recorded noise, integration before absolute-value
% pooling, and postcomparison Gaussian noise are retained. No lapse or bias
% parameter is included by default. Optional LapseMode='session' profiles a
% random-guess lapse fraction separately per session; 'fixed' specifies one
% fraction for all sessions. P(correct)=lambda/2+(1-lambda)*Phi(D/sigma).
% Thus lambda=.10 means 10% guessing trials and a 95% performance ceiling,
% not a 10% asymptotic error rate. No choice-bias parameter is fitted.
% JointAll is a fitted score, NOT held-out accuracy. Lapse estimates are
% model-dependent and can be poorly identified without easy trials.
%
% Temporal kernels still use other sessions' 100:0 choices only. Optional
% Session100 fits sigma on same-session 100:0 trials and scores divided
% trials. Do not compare the absolute LLs of the two schemes: their scoring
% sets differ. The old OtherAll/common-training-sigma scheme is removed.
%
% BOOTSTRAP: resample whole sessions, refit the leave-one-session-out kernels
% and all sigma/lapse/q profiles, and take percentile limits on q. Every copy of
% the evaluated source session is excluded from its kernel training. Draws
% containing only one distinct session, or without an estimable kernel, are
% rejected and counted. The nominal 95% interval is on the supplied q grid;
% coverage is not guaranteed, especially at q=0/1 and with few sessions.
% bootstrapFractionAtMax is a resampling frequency, NOT P(max is true).
%
% RECOVERY (optional and expensive): generate new Bernoulli choices on the
% recorded stimuli from specified true q values, using the original fitted
% per-session kernels and sigma at each generating q. Refit the complete
% pipeline; optionally bootstrap each recovery replicate to check coverage.
% Recovery uses the selected lapse mode both to generate and fit, conditional
% on these stimuli and generating kernels/sigmas/lapses. It does not test
% neural mechanism uniqueness or misspecification of the lapse process.
%
% Small recovery smoke run:
% R = pilotIDLPNormProfile('Animal','Neesha','NBootstrap',100, ...
%     'NRecovery',10,'RecoveryQ',[0 .5 1],'RecoveryNBootstrap',50);
% Increase NRecovery/RecoveryNBootstrap for a substantive coverage check.
% RecoveryNBootstrap=0 checks bias/endpoint frequency only (coverage is NaN).
%
% Main outputs: summaryTable, profileLL, sessionLLByQScheme,
% sigmaBySessionQScheme, lapseRateBySessionQScheme, bootstrapQHat, bootstrapProfileLL,
% allocationFitTable (observed/predicted accuracy at the fitted primary q),
% and recovery.summaryTable. Results/PDF use animal-specific filenames; lapse
% runs add _LapseSession or _LapseFixedXXpc so the no-lapse run is preserved.
% Requires the existing domainFolder/selectAnalysisFiles/validFolder helpers.

p=inputParser;
p.KeepUnmatched=true;
p.PartialMatching=false;
p.addParameter('TwoAllocationMode',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('Visible','on',@(x) any(strcmpi(string(x),["on" "off"])));
p.addParameter('SessionFolder','',@(x) ischar(x)||isstring(x));
p.addParameter('Animal','Meetz',@(x) (ischar(x)||isstring(x))&& ...
  isscalar(string(x))&&strlength(string(x))>0);
p.addParameter('QGrid',0:.025:1,@(x) isnumeric(x)&&isvector(x)&& ...
  numel(x)>=3&&all(isfinite(x))&&all(x>=0)&&all(x<=1));
p.addParameter('NBootstrap',1000,@isNonnegativeInteger);
p.addParameter('IncludeSession100',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('LapseMode','none',@(x) (ischar(x)||isstring(x))&& ...
  isscalar(string(x))&&any(strcmpi(string(x),["none" "fixed" "session"])));
p.addParameter('LapseRate',0,@isLapseFraction);
p.addParameter('MaxLapseRate',.4,@isLapseFraction);
p.addParameter('Seed',14127,@(x) isnumeric(x)&&isscalar(x)&& ...
  isfinite(x)&&x>=0&&x==round(x)&&x<=2^32-1);
p.addParameter('NRecovery',0,@isNonnegativeInteger);
p.addParameter('RecoveryQ',[0 .5 1],@(x) isnumeric(x)&&isvector(x)&& ...
  ~isempty(x)&&all(isfinite(x))&&all(x>=0)&&all(x<=1));
p.addParameter('RecoveryNBootstrap',0,@isNonnegativeInteger);
p.addParameter('Verbose',true,@(x) islogical(x)&&isscalar(x));
p.addParameter('SaveResult',true,@(x) islogical(x)&&isscalar(x));
p.addParameter('ExportPDF',true,@(x) islogical(x)&&isscalar(x));
p.parse(varargin{:});
opts=p.Results;
opts.LapseMode=lower(char(string(opts.LapseMode)));
animal=char(string(opts.Animal));
q=unique(double(opts.QGrid(:))');
assert(all(ismember([0 .5 1],q)), ...
  'pilotIDLPNormProfile:RequiredQ', ...
  'QGrid must contain 0, 0.5, and 1 for the prespecified comparisons.');
root=domainFolder(mfilename('fullpath'));
folder=char(string(opts.SessionFolder));
if isempty(folder), folder=fullfile(root,'Data','SessionAnalysis'); end
assert(isfolder(folder),'pilotIDLPNormProfile:MissingFolder', ...
  'SessionAnalysis folder not found: %s',folder);
selection=p.Unmatched;
fields=fieldnames(selection);
args=cell(1,2*numel(fields));
for k=1:numel(fields)
  args{2*k-1}=fields{k}; args{2*k}=selection.(fields{k});
end
paths=selectAnalysisFiles(folder,'Animal',animal, ...
  'TwoAllocationMode',opts.TwoAllocationMode,args{:});
if ischar(paths), paths={paths}; end
if isstring(paths), paths=cellstr(paths(:)); end
paths=paths(:);
nS=numel(paths);
assert(nS>=3,'pilotIDLPNormProfile:TooFewSessions', ...
  'At least three distinct sessions are required.');
A=cell(nS,1);
ids=strings(nS,1);
ratios=nan(nS,1);
for s=1:nS
  loaded=load(paths{s},'sessionAnalysis');
  assert(isfield(loaded,'sessionAnalysis'), ...
    'pilotIDLPNormProfile:MissingAnalysis','%s lacks sessionAnalysis.',paths{s});
  idlAssertAnalysisMode(loaded.sessionAnalysis,opts.TwoAllocationMode);
  validateSession(loaded.sessionAnalysis,animal);
  A{s}=prepareSession(loaded.sessionAnalysis);
  ids(s)=string(loaded.sessionAnalysis.sourceSessionID);
  ratios(s)=loaded.sessionAnalysis.sessionHeader.splitUnsplitRatio;
end
assert(numel(unique(ids))==nS,'pilotIDLPNormProfile:DuplicateSessions', ...
  'Selection contains duplicate sourceSessionID values.');
times=A{1}.times;
assert(all(cellfun(@(x) isequal(x.times,times),A)), ...
  'pilotIDLPNormProfile:TimeGrid', ...
  'Sessions have incompatible video-frame time grids.');
schemeName="JointAll";
if opts.IncludeSession100, schemeName=[schemeName;"Session100"]; end
nSchemes=numel(schemeName);
counts=ones(nS,1);
summaries=kernelSummaries(A);
F=fitProfiles(A,summaries,counts,q,opts.IncludeSession100,opts);
assert(F.ok,'pilotIDLPNormProfile:KernelOutcomes', ...
  'Other-session 100:0 trials must provide a nonzero correct/error kernel.');
profile=F.profileLL;
relativeProfile=profile-max(profile,[],1);
[~,best]=max(profile,[],1);
qHat=q(best)';

% Restore the caller's RNG state even on failure. Recovery uses a separate
% stream so changing NBootstrap does not change simulated recovery choices.
oldRng=rng;
cleanup=onCleanup(@() rng(oldRng)); %#ok<NASGU>
rng(opts.Seed,'twister');
B=bootstrapProfiles(A,summaries,q,opts.IncludeSession100, ...
  opts.NBootstrap,opts.Verbose,opts);
[low,high,massAtMax,massAtLinear]=bootstrapLimits(B.qHat,nSchemes);

pHat=reciprocalQ(qHat);
lower95Q=low'; upper95Q=high';
lower95P=reciprocalQ(upper95Q); upper95P=reciprocalQ(lower95Q);
relativeLLatMax=relativeProfile(q==0,:)';
relativeLLatL2=relativeProfile(q==.5,:)';
relativeLLatL1=relativeProfile(q==1,:)';
bootstrapFractionAtMax=massAtMax';
bootstrapFractionAtLinear=massAtLinear';
meanSessionLapseRate=nan(nSchemes,1);
maxSessionLapseRate=nan(nSchemes,1);
for k=1:nSchemes
  rates=F.lapseRate(:,best(k),k);
  meanSessionLapseRate(k)=mean(rates);
  maxSessionLapseRate(k)=max(rates);
end
summaryTable=table(schemeName,qHat,pHat,lower95Q,upper95Q, ...
  lower95P,upper95P,relativeLLatMax,relativeLLatL2,relativeLLatL1, ...
  bootstrapFractionAtMax,bootstrapFractionAtLinear, ...
  meanSessionLapseRate,maxSessionLapseRate);
result=struct();
result.analysisVersion=4;
result.twoAllocationMode=opts.TwoAllocationMode;
result.splitUnsplitRatioBySession=ratios;
result.createdAt=datetime('now');
result.animal=animal;
result.options=opts;
result.qGrid=q;
result.pGrid=reciprocalQ(q);
result.schemeNames=schemeName;
result.sourceFiles=string(paths);
result.sourceSessionID=ids;
result.nTrials=cellfun(@(x) numel(x.y),A);
result.nDividedTrials=cellfun(@(x) sum(~x.is100),A);
result.nScoredTrialsBySessionScheme=result.nTrials;
if opts.IncludeSession100
  result.nScoredTrialsBySessionScheme=[result.nTrials result.nDividedTrials];
end
result.stepTimeMS=times;
result.looKernel=F.kernels;
result.sigmaBySessionQScheme=F.sigma;
result.lapseRateBySessionQScheme=F.lapseRate;
result.sessionLLByQScheme=F.sessionLL;
% Compatibility: this contains ONLY the actual Session100 divided scores.
% It is empty by default; JointAll LLs must not masquerade as held-out LLs.
result.heldOutDividedLL=[];
if opts.IncludeSession100, result.heldOutDividedLL=F.sessionLL(:,:,2); end
result.profileLL=profile;
result.relativeProfileLL=relativeProfile;
result.bootstrapQHat=B.qHat;
result.bootstrapProfileLL=B.profileLL;
result.bootstrapSessionCounts=B.counts;
result.bootstrapRejectedDraws=B.rejectedDraws;
result.bootstrapSeed=opts.Seed;
result.summaryTable=summaryTable;
result.allocationFitTable=allocationFit(A,F,q,best(1));
result.modelContract=['Absolute coherence including pedestal and recorded ' ...
  'noise -> temporal integration within stream -> absolute-value p-norm ' ...
  'within patch -> Change minus No-change -> Phi(D/sigma); optional ' ...
  'random-guess mixture lambda/2+(1-lambda)*Phi(D/sigma). No choice bias.'];
result.bootstrapDescription=['Whole-session percentile bootstrap; refit ' ...
  'LOO kernels and per-session sigma/lapse/q profiles; exclude all copies of ' ...
  'evaluated source session. Nominal 95% limits on q grid; coverage must ' ...
  'be checked, especially at boundaries. Not a posterior model probability.'];
result.schemeDescription=['JointAll: one pooled q, sigma profiled separately ' ...
  'per session using all allocations; fitted LL. Optional Session100: ' ...
  'same-session 100:0 sigma, divided-trial prediction score.'];
result.recovery=struct();
if opts.NRecovery>0
  recoverySeed=mod(double(opts.Seed)+104729,2^32);
  rng(recoverySeed,'twister');
  result.recovery=runRecovery(A,F,q,opts,recoverySeed);
end
if opts.SaveResult
  dest=validFolder(fullfile(root,'Data','AcrossSessions'));
  save(fullfile(dest,['PNormPilot_',animal,idlAllocationTag(opts.TwoAllocationMode),lapseFileTag(opts),'.mat']),'result','-v7.3');
end
fig=plotPilot(result);
if opts.ExportPDF
  dest=validFolder(fullfile(root,'Plots','AcrossSessions'));
  exportgraphics(fig,fullfile(dest,['PNormPilot_',animal,idlAllocationTag(opts.TwoAllocationMode),lapseFileTag(opts),'.pdf']), ...
    'ContentType','vector');
end
disp(summaryTable)
if opts.NRecovery>0, disp(result.recovery.summaryTable), end
end

function tf=isNonnegativeInteger(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==round(x);
end

function tf=isLapseFraction(x)
tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<1;
end

function tag=lapseFileTag(opts)
tag='';
if strcmp(opts.LapseMode,'session')
  tag='_LapseSession';
elseif strcmp(opts.LapseMode,'fixed')
  tag=['_LapseFixed',strrep(sprintf('%g',100*opts.LapseRate),'.','p'),'pc'];
end
end

function validateSession(A,animal)
required={'sessionHeader','sourceSessionID','trialTable','stepTimeMS', ...
  'changeNoiseByFrameTrial','noChangeNoiseByFrameTrial'};
assert(all(isfield(A,required)), ...
  'pilotIDLPNormProfile:MissingFields','Session analysis lacks required fields.');
assert(isfield(A.sessionHeader,'animal')&& ...
  strcmpi(string(A.sessionHeader.animal),string(animal)), ...
  'pilotIDLPNormProfile:Animal','Selected session has a different animal.');
T=A.trialTable;
names={'correct','alloc0','alloc1','totalStepCohPC','prestepCohPC', ...
  'unsignedAllocPC','highStreamIndex'};
assert(istable(T)&&all(ismember(names,T.Properties.VariableNames)), ...
  'pilotIDLPNormProfile:TrialFields','Required trialTable fields are missing.');
n=height(T); f=numel(A.stepTimeMS);
assert(n>0&&f>0&&all(isfinite(double(A.stepTimeMS(:))))&& ...
  isequal(size(A.changeNoiseByFrameTrial),[f n 2])&& ...
  isequal(size(A.noChangeNoiseByFrameTrial),[f n 2]), ...
  'pilotIDLPNormProfile:NoiseShape','Unexpected time grid or noise-array dimensions.');
% The high stream is undefined for ties; only 100:0 indices are needed.
for k=1:numel(names)-1
  v=double(T.(names{k})(:));
  assert(numel(v)==n&&all(isfinite(v)), ...
    'pilotIDLPNormProfile:TrialValues','Invalid trial field %s.',names{k});
end
assert(all(ismember(double(T.correct(:)),[0 1]))&& ...
  all(double(T.alloc0(:))>=0)&&all(double(T.alloc1(:))>=0)&& ...
  all(abs(double(T.alloc0(:))+double(T.alloc1(:))-1)<1e-6)&& ...
  all(double(T.totalStepCohPC(:))>=0)&&all(double(T.prestepCohPC(:))>=0), ...
  'pilotIDLPNormProfile:TrialValues','Invalid choices, allocations, or coherences.');
assert(all(isfinite(A.changeNoiseByFrameTrial(:)))&& ...
  all(isfinite(A.noChangeNoiseByFrameTrial(:))), ...
  'pilotIDLPNormProfile:NoiseValues','Noise arrays contain nonfinite values.');
is100=abs(double(T.unsignedAllocPC(:))-100)<1e-6;
assert(numel(T.highStreamIndex)==n&&any(is100)&& ...
  all(ismember(double(T.highStreamIndex(is100)),[1 2])), ...
  'pilotIDLPNormProfile:ReferenceStream','100:0 reference trials/streams are invalid.');
end

function S=prepareSession(A)
T=A.trialTable;
S=struct();
S.times=double(A.stepTimeMS(:));
S.y=double(T.correct(:));
S.is100=abs(double(T.unsignedAllocPC(:))-100)<1e-6;
S.alloc=double(T.unsignedAllocPC(:));
S.step=double(T.totalStepCohPC(:));
S.signal=[double(T.alloc0(:)) double(T.alloc1(:))].*S.step;
S.pedestal=double(T.prestepCohPC(:));
S.Cnoise=double(A.changeNoiseByFrameTrial);
S.Nnoise=double(A.noChangeNoiseByFrameTrial);
rows=find(S.is100);
S.calNoise=zeros(numel(S.times),numel(rows));
for j=1:numel(rows)
  stream=double(T.highStreamIndex(rows(j)));
  S.calNoise(:,j)=S.Cnoise(:,rows(j),stream);
end
end

function H=kernelSummaries(A)
nS=numel(A); f=numel(A{1}.times);
H=struct('correctSum',zeros(f,nS),'errorSum',zeros(f,nS), ...
  'nCorrect',zeros(nS,1),'nError',zeros(nS,1));
H.hasNoise=cellfun(@(x) any(abs(x.calNoise)>1e-8,'all'),A);
for s=1:nS
  good=logical(A{s}.y(A{s}.is100));
  H.correctSum(:,s)=sum(A{s}.calNoise(:,good),2);
  H.errorSum(:,s)=sum(A{s}.calNoise(:,~good),2);
  H.nCorrect(s)=sum(good);
  H.nError(s)=sum(~good);
end
end

function [w,ok]=trainingKernel(H,counts,held)
% Multiplicities weight training statistics. Remove EVERY copy of held.
trainCounts=counts;
trainCounts(held)=0;
if ~any(H.hasNoise(:) & trainCounts(:)>0)
  % Constant noise-free steps need no estimated temporal weighting.
  w=ones(size(H.correctSum,1),1)/size(H.correctSum,1); ok=true; return
end
nCorrect=H.nCorrect'*trainCounts; nError=H.nError'*trainCounts;
w=nan(size(H.correctSum,1),1);
ok=nCorrect>0&&nError>0;
if ~ok, return, end
w=(H.correctSum*trainCounts)/nCorrect-(H.errorSum*trainCounts)/nError;
if sum(w)<0, w=-w; end
mass=sum(abs(w));
ok=isfinite(mass)&&mass>0;
if ok, w=w/mass; end
end

function F=fitProfiles(A,H,counts,q,includeSession100,opts)
nS=numel(A); nQ=numel(q); nSchemes=1+double(includeSession100);
F=struct('ok',true,'sessionLL',nan(nS,nQ,nSchemes), ...
  'sigma',nan(nS,nQ,nSchemes),'lapseRate',nan(nS,nQ,nSchemes), ...
  'kernels',nan(numel(A{1}.times),nS), ...
  'profileLL',zeros(nQ,nSchemes));
for s=find(counts(:)>0)'
  [w,ok]=trainingKernel(H,counts,s);
  if ~ok, F.ok=false; return, end
  F.kernels(:,s)=w;
  [C,N]=integratedStreams(A{s},w);
  y=A{s}.y;
  for j=1:nQ
    D=evidence(C,N,q(j));
    [sigma,lambda]=fitDecisionParameters(D,y,opts);
    F.sigma(s,j,1)=sigma;
    F.lapseRate(s,j,1)=lambda;
    F.sessionLL(s,j,1)=bernoulliLL(lapseProbability(D,sigma,lambda),y);
    if includeSession100
      cal=A{s}.is100; test=~cal;
      [sigma,lambda]=fitDecisionParameters(D(cal),y(cal),opts);
      F.sigma(s,j,2)=sigma;
      F.lapseRate(s,j,2)=lambda;
      F.sessionLL(s,j,2)=bernoulliLL(lapseProbability(D(test),sigma,lambda),y(test));
    end
  end
  F.profileLL=F.profileLL+counts(s)*reshape(F.sessionLL(s,:,:),nQ,nSchemes);
end
end

function [C,N]=integratedStreams(S,w)
n=numel(S.y);
C=zeros(n,2); N=zeros(n,2);
for stream=1:2
  % Algebraically identical to adding pedestal/step at every video frame.
  C(:,stream)=S.Cnoise(:,:,stream)'*w+(S.pedestal+S.signal(:,stream))*sum(w);
  N(:,stream)=S.Nnoise(:,:,stream)'*w+S.pedestal*sum(w);
end
end

function D=evidence(C,N,q)
D=poolStreams(C,q)-poolStreams(N,q);
end

function v=poolStreams(X,q)
X=abs(X);
m=max(X,[],2);
if q==0, v=m; return, end
% Scale before exponentiation: tiny positive q must not overflow, and an
% all-zero row must remain zero rather than produce 0/0.
denom=m; denom(denom==0)=1;
v=m.*sum((X./denom).^(1/q),2).^q;
end

function sigma=fitSigma(D,y,lambda)
if nargin<3, lambda=0; end
scale=max(std(D),1);
objective=@(z) -bernoulliLL(lapseProbability(D,exp(z),lambda),y);
bounds=log([scale*1e-3 scale*1e3]);
z=fminbnd(objective,bounds(1),bounds(2), ...
  optimset('Display','off','TolX',1e-8));
sigma=exp(z);
end

function [sigma,lambda]=fitDecisionParameters(D,y,opts)
lambda=0;
if strcmp(opts.LapseMode,'fixed'), lambda=opts.LapseRate; end
sigma=fitSigma(D,y,lambda);
if ~strcmp(opts.LapseMode,'session')||opts.MaxLapseRate==0, return, end
% Explicit endpoint candidates prevent logistic parametrization from
% excluding zero lapses or the configured upper bound. Interior solutions
% use multiple starts; sigma is constrained to the existing broad bounds.
bestLL=bernoulliLL(lapseProbability(D,sigma,lambda),y);
cap=opts.MaxLapseRate;
capSigma=fitSigma(D,y,cap);
capLL=bernoulliLL(lapseProbability(D,capSigma,cap),y);
if capLL>bestLL+1e-7, sigma=capSigma; lambda=cap; bestLL=capLL; end
scale=max(std(D),1); bounds=log([scale*1e-3 scale*1e3]);
for guess=cap*[.05 .25 .625]
  startSigma=fitSigma(D,y,guess);
  start=[log(startSigma),log(guess/(cap-guess))];
  objective=@(v) lapseObjective(v,D,y,cap,bounds);
  [v,nll]=fminsearch(objective,start,optimset('Display','off', ...
    'TolX',1e-7,'TolFun',1e-7,'MaxIter',1500,'MaxFunEvals',3000));
  if isfinite(nll)&&-nll>bestLL+1e-7
    sigma=exp(v(1)); lambda=boundedLapse(v(2),cap); bestLL=-nll;
  end
end
end

function nll=lapseObjective(v,D,y,cap,bounds)
if any(~isfinite(v))||v(1)<bounds(1)||v(1)>bounds(2)
  nll=Inf; return
end
nll=-bernoulliLL(lapseProbability(D,exp(v(1)),boundedLapse(v(2),cap)),y);
end

function lambda=boundedLapse(z,cap)
if z>=0
  lambda=cap/(1+exp(-z));
else
  ez=exp(z); lambda=cap*ez/(1+ez);
end
end

function p=lapseProbability(D,sigma,lambda)
p=lambda/2+(1-lambda)*normalCDF(D/sigma);
p=min(max(p,1e-12),1-1e-12);
end

function p=normalCDF(D)
p=.5*erfc(-double(D)/sqrt(2));
p=min(max(p,1e-12),1-1e-12);
end

function ll=bernoulliLL(p,y)
p=min(max(double(p(:)),1e-12),1-1e-12);
y=double(y(:));
ll=sum(y.*log(p)+(1-y).*log(1-p));
end

function B=bootstrapProfiles(A,H,q,includeSession100,nBoot,verbose,opts)
nS=numel(A); nQ=numel(q); nSchemes=1+double(includeSession100);
B=struct('qHat',nan(nBoot,nSchemes),'profileLL',nan(nBoot,nQ,nSchemes), ...
  'counts',zeros(nBoot,nS),'rejectedDraws',0);
for b=1:nBoot
  accepted=false;
  for attempt=1:10000
    counts=accumarray(randi(nS,nS,1),1,[nS 1]);
    if sum(counts>0)<2
      B.rejectedDraws=B.rejectedDraws+1; continue
    end
    F=fitProfiles(A,H,counts,q,includeSession100,opts);
    if ~F.ok
      B.rejectedDraws=B.rejectedDraws+1; continue
    end
    accepted=true;
    [~,best]=max(F.profileLL,[],1);
    B.qHat(b,:)=q(best);
    B.profileLL(b,:,:)=reshape(F.profileLL,1,nQ,nSchemes);
    B.counts(b,:)=counts';
    break
  end
  assert(accepted,'pilotIDLPNormProfile:BootstrapKernel', ...
    'Cannot obtain an estimable bootstrap kernel after 10000 draws.');
  if verbose&&(b==1||mod(b,25)==0||b==nBoot)
    fprintf('Full-pipeline bootstrap %d/%d (%d rejected draws)\n', ...
      b,nBoot,B.rejectedDraws);
    drawnow;
  end
end
end

function [low,high,massAtMax,massAtLinear]=bootstrapLimits(bootQ,nSchemes)
low=nan(1,nSchemes); high=low; massAtMax=low; massAtLinear=low;
n=size(bootQ,1);
if n==0, return, end
for k=1:nSchemes
  sorted=sort(bootQ(:,k));
  low(k)=sorted(max(1,ceil(.025*n)));
  high(k)=sorted(min(n,ceil(.975*n)));
  massAtMax(k)=mean(bootQ(:,k)==0);
  massAtLinear(k)=mean(bootQ(:,k)==1);
end
end

function p=reciprocalQ(q)
p=1./q;
p(q==0)=Inf;
end

function T=allocationFit(A,F,q,best)
% Descriptive check using each allocation's actual sampled step grid. These
% are fitted accuracies, not thresholds or out-of-sample predictions.
alloc=[]; step=[]; y=[]; probability=[];
for s=1:numel(A)
  [C,N]=integratedStreams(A{s},F.kernels(:,s));
  D=evidence(C,N,q(best));
  alloc=[alloc;A{s}.alloc]; %#ok<AGROW>
  step=[step;A{s}.step]; %#ok<AGROW>
  y=[y;A{s}.y]; %#ok<AGROW>
  probability=[probability;lapseProbability(D,F.sigma(s,best,1), ...
    F.lapseRate(s,best,1))]; %#ok<AGROW>
end
levels=unique(alloc);
nTrials=zeros(size(levels)); nCorrect=nTrials; meanStepCohPC=nTrials;
observedPCorrect=nTrials; fittedPCorrect=nTrials;
for k=1:numel(levels)
  rows=alloc==levels(k);
  nTrials(k)=sum(rows); nCorrect(k)=sum(y(rows));
  meanStepCohPC(k)=mean(step(rows));
  observedPCorrect(k)=mean(y(rows)); fittedPCorrect(k)=mean(probability(rows));
end
unsignedAllocPC=levels;
T=table(unsignedAllocPC,nTrials,nCorrect,meanStepCohPC, ...
  observedPCorrect,fittedPCorrect);
end

function R=runRecovery(A,F,q,opts,seed)
trueQ=unique(double(opts.RecoveryQ(:))');
nTruth=numel(trueQ); nR=opts.NRecovery; nS=numel(A);
R=struct();
R.seed=seed;
R.trueQ=trueQ;
R.qHat=nan(nR,nTruth);
R.lower95Q=nan(nR,nTruth); R.upper95Q=nan(nR,nTruth);
R.generatingSigma=nan(nS,nTruth);
R.generatingLapseRate=nan(nS,nTruth);
R.bootstrapRejectedDraws=zeros(nR,nTruth);
R.nBootstrapPerReplicate=opts.RecoveryNBootstrap;
R.description=['Conditional parametric recovery on recorded stimuli; ' ...
  'generate from original LOO kernels and per-session sigma/lapse, re-estimate ' ...
  '100:0 kernels and JointAll sigma/lapse/q using selected lapse mode; optional full-pipeline session ' ...
  'bootstrap per replicate. Coverage estimates require many replicates.'];
generationProbability=cell(nS,nTruth);
for s=1:nS
  [C,N]=integratedStreams(A{s},F.kernels(:,s));
  for t=1:nTruth
    D=evidence(C,N,trueQ(t));
    [sigma,lambda]=fitDecisionParameters(D,A{s}.y,opts);
    R.generatingSigma(s,t)=sigma;
    R.generatingLapseRate(s,t)=lambda;
    generationProbability{s,t}=lapseProbability(D,sigma,lambda);
  end
end
% Separate choice/interval RNG states. Increasing RecoveryNBootstrap leaves
% the simulated datasets unchanged, allowing fair numerical comparisons.
choiceState=rng;
rng(mod(double(seed)+13007,2^32),'twister');
intervalState=rng;
for t=1:nTruth
  for r=1:nR
    rng(choiceState);
    simulated=A;
    for s=1:nS
      simulated{s}.y=double(rand(size(A{s}.y))<generationProbability{s,t});
    end
    choiceState=rng;
    H=kernelSummaries(simulated);
    fitted=fitProfiles(simulated,H,ones(nS,1),q,false,opts);
    assert(fitted.ok,'pilotIDLPNormProfile:RecoveryKernel', ...
      'Recovery dataset has no estimable LOO kernel; inspect sampling.');
    [~,idx]=max(fitted.profileLL);
    R.qHat(r,t)=q(idx);
    if opts.RecoveryNBootstrap>0
      rng(intervalState);
      B=bootstrapProfiles(simulated,H,q,false,opts.RecoveryNBootstrap,false,opts);
      intervalState=rng;
      [lo,hi]=bootstrapLimits(B.qHat,1);
      R.lower95Q(r,t)=lo; R.upper95Q(r,t)=hi;
      R.bootstrapRejectedDraws(r,t)=B.rejectedDraws;
    end
    if opts.Verbose
      fprintf('Recovery true q=%.3f: %d/%d; fitted q=%.3f\n', ...
        trueQ(t),r,nR,R.qHat(r,t));
      drawnow;
    end
  end
end
qTrue=trueQ'; pTrue=reciprocalQ(qTrue);
nReplicates=repmat(nR,nTruth,1);
meanRecoveredQ=mean(R.qHat,1)';
biasQ=meanRecoveredQ-qTrue;
rmseQ=sqrt(mean((R.qHat-trueQ).^2,1))';
fractionAtMax=mean(R.qHat==0,1)';
fractionAtLinear=mean(R.qHat==1,1)';
coverage95=nan(nTruth,1); coverageMonteCarloSE=coverage95;
meanIntervalWidthQ=coverage95;
if opts.RecoveryNBootstrap>0
  covered=R.lower95Q<=trueQ & R.upper95Q>=trueQ;
  coverage95=mean(covered,1)';
  coverageMonteCarloSE=sqrt(coverage95.*(1-coverage95)/nR);
  meanIntervalWidthQ=mean(R.upper95Q-R.lower95Q,1)';
end
R.summaryTable=table(qTrue,pTrue,nReplicates,meanRecoveredQ,biasQ,rmseQ, ...
  fractionAtMax,fractionAtLinear,coverage95,coverageMonteCarloSE,meanIntervalWidthQ);
end

function fig=plotPilot(R)
fig=figure('Name',sprintf('IDL p-norm pilot | %s',R.animal), ...
  'NumberTitle','off','Color','w','Position',[100 100 1350 720], ...
  'Visible',char(R.options.Visible));
if strcmpi(R.options.Visible,'on'), fig.WindowStyle='docked'; end
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
sgtitle(tl,['IDL p-norm | ',char(R.animal),idlAllocationTag(R.twoAllocationMode)], ...
  'Interpreter','none');
colors=[.12 .42 .78; .83 .35 .17];
nSchemes=numel(R.schemeNames);
for scheme=1:nSchemes
  ax=nexttile(tl,scheme); hold(ax,'on');
  y=R.relativeProfileLL(:,scheme);
  plot(ax,R.qGrid,y,'-','Color',colors(scheme,:),'LineWidth',2);
  [~,idx]=max(y);
  plot(ax,R.qGrid(idx),y(idx),'o','MarkerFaceColor',colors(scheme,:), ...
    'MarkerEdgeColor',colors(scheme,:),'MarkerSize',7);
  xline(ax,0,':','Max'); xline(ax,.5,':','L2'); xline(ax,1,':','L1');
  yline(ax,0,':','HandleVisibility','off');
  xlabel(ax,'q = 1/p');
  if scheme==1
    ylabel(ax,'All-allocation fitted LL - peak LL');
  else
    ylabel(ax,'Divided-trial prediction LL - peak LL');
  end
  if isempty(R.bootstrapQHat)
    title(ax,sprintf('%s: peak q = %.3f',R.schemeNames(scheme), ...
      R.summaryTable.qHat(scheme)));
  else
    xline(ax,R.summaryTable.lower95Q(scheme),'--','HandleVisibility','off');
    xline(ax,R.summaryTable.upper95Q(scheme),'--','HandleVisibility','off');
    title(ax,sprintf('%s: q = %.3f; bootstrap interval [%.3f, %.3f]', ...
      R.schemeNames(scheme),R.summaryTable.qHat(scheme), ...
      R.summaryTable.lower95Q(scheme),R.summaryTable.upper95Q(scheme)));
  end
  xlim(ax,[0 1]); grid(ax,'on'); box(ax,'off');
  ax=nexttile(tl,scheme+2); hold(ax,'on');
  d=R.sessionLLByQScheme(:,:,scheme)- ...
    R.sessionLLByQScheme(:,R.qGrid==.5,scheme);
  % imagesc ignores interior x coordinates. Use grid indices and truthful
  % tick labels for a nonuniform user-supplied grid.
  uniformGrid=max(abs(diff(R.qGrid)-mean(diff(R.qGrid))))<1e-10;
  if uniformGrid
    imagesc(ax,R.qGrid,1:numel(R.sourceSessionID),d);
    xline(ax,.5,'w:','LineWidth',1.5); xlim(ax,[0 1]);
    xlabel(ax,'q = 1/p');
  else
    imagesc(ax,1:numel(R.qGrid),1:numel(R.sourceSessionID),d);
    ticks=unique([1:ceil(numel(R.qGrid)/6):numel(R.qGrid), ...
      find(R.qGrid==.5),numel(R.qGrid)]);
    xticks(ax,ticks); xticklabels(ax,arrayfun(@(q) sprintf('%.3g',q),R.qGrid(ticks),'UniformOutput',false));
    xline(ax,find(R.qGrid==.5),'w:','LineWidth',1.5);
    xlim(ax,[.5 numel(R.qGrid)+.5]); xlabel(ax,'q = 1/p (grid points)');
  end
  colormap(ax,parula); colorbar(ax);
  ylim(ax,[.5 numel(R.sourceSessionID)+.5]); ylabel(ax,'Session');
  title(ax,sprintf('%s: session LL relative to L2',R.schemeNames(scheme)));
end
if nSchemes==1
  ax=nexttile(tl,2); hold(ax,'on');
  if isempty(R.bootstrapQHat)
    text(ax,.5,.5,'Bootstrap disabled','HorizontalAlignment','center');
  else
    % Bar plot handles custom/nonuniform grids without misleading bins.
    frequency=arrayfun(@(x) mean(R.bootstrapQHat(:,1)==x),R.qGrid);
    bar(ax,R.qGrid,frequency,'FaceColor',colors(1,:),'EdgeColor','none');
    xline(ax,R.summaryTable.qHat(1),'k--');
  end
  xlim(ax,[-.025 1.025]); xlabel(ax,'Bootstrap best q');
  ylabel(ax,'Resampling fraction'); title(ax,'Full-pipeline session bootstrap');
  grid(ax,'on'); box(ax,'off');
  ax=nexttile(tl,4);
  [~,idx]=max(R.profileLL(:,1));
  bar(ax,R.sigmaBySessionQScheme(:,idx,1),'FaceColor',colors(1,:));
  xlabel(ax,'Session'); ylabel(ax,'Decision sigma (coherence percentage points)');
  title(ax,'Session sigma at pooled best q'); grid(ax,'on'); box(ax,'off');
  if ~strcmp(R.options.LapseMode,'none')
    yyaxis(ax,'right');
    plot(ax,100*R.lapseRateBySessionQScheme(:,idx,1),'o-', ...
      'Color',colors(2,:),'LineWidth',1.5);
    ylabel(ax,'Guessing trials (%)');
    title(ax,'Session sigma and lapse fraction at pooled best q');
  end
end
sgtitle(tl,sprintf(['%s | p-norm profile | %d sessions, %d total trials | ' ...
  'lapses: %s | %d bootstrap replicates'],R.animal,numel(R.sourceSessionID), ...
  sum(R.nTrials),R.options.LapseMode,size(R.bootstrapQHat,1)));
end
