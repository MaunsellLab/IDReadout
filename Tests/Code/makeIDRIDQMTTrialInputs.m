function inputs = makeIDRIDQMTTrialInputs( ...
  idqSummaryPath,idrSummaryPath,idqPedestalPC)
% makeIDRIDQMTTrialInputs  Standardize actual IDQ/IDR step-period stimuli.
%
% Output conditions use drift-relative directions. componentCoherence is
% the actual step-period coherence distribution (percent coherence), while
% signalCoherence is the corresponding noise-free distribution.

arguments
  idqSummaryPath (1,1) string
  idrSummaryPath (1,1) string
  idqPedestalPC double {mustBeFinite,mustBeNonnegative}
end

Q = load(idqSummaryPath,'acrossSummary');
R = load(idrSummaryPath,'inventory');
if ~isfield(Q,'acrossSummary') || ~isfield(R,'inventory')
  error('makeIDRIDQMTTrialInputs:MissingSummary', ...
    'Expected acrossSummary and inventory in the IDQ and IDR files.');
end
inputs = [makeIDQ(Q.acrossSummary,idqPedestalPC); makeIDR(R.inventory)];
end

function condition = makeIDQ(S,idqPedestalPC)
if ~isfield(S,'trialTable')
  error('makeIDRIDQMTTrialInputs:MissingIDQTrialTable', ...
    'IDQ acrossSummary has no trialTable.');
end
T = S.trialTable;
vars = { ...
  'sessionIndex','trialIndex','dirIndex','stepCoh','hasStepNoise','correct', ...
  'rectStepChangeNoisePredDir1','rectStepChangeNoisePredDir2', ...
  'rectStepChangeNoisePredDir3','rectStepNoChangeNoisePredDir1', ...
  'rectStepNoChangeNoisePredDir2','rectStepNoChangeNoisePredDir3'};
requireTableVariables(T,vars,'IDQ');
use = logical(T.hasStepNoise) & isfinite(T.stepCoh) & ...
  ismember(double(T.dirIndex),1:3);
T = T(use,:);

changeAbs = tableColumns(T,vars(7:9));
noChangeAbs = tableColumns(T,vars(10:12));
finiteRows = all(isfinite(changeAbs),2) & all(isfinite(noChangeAbs),2);
nExcluded = sum(~finiteRows);
T = T(finiteRows,:); changeAbs=changeAbs(finiteRows,:);
noChangeAbs=noChangeAbs(finiteRows,:);
n = height(T);
pedestalByTrial = resolveIDQPedestal(idqPedestalPC,T.sessionIndex);
changeNoise = nan(n,3); noChangeNoise = nan(n,3);
for i = 1:n
  drift = double(T.dirIndex(i));
  relativeOrder = [drift, mod(drift,3)+1, mod(drift-2,3)+1];
  changeNoise(i,:) = changeAbs(i,relativeOrder);
  noChangeNoise(i,:) = noChangeAbs(i,relativeOrder);
end
step = double(T.stepCoh);
noChangeSignal = repmat(pedestalByTrial,1,3);
changeSignal = noChangeSignal;
changeSignal(:,1) = changeSignal(:,1)+step;

condition = baseCondition("IDQ",120,[0 120 -120], ...
  ["primary","probePlus","probeMinus"],T.sessionIndex,T.trialIndex);
condition.changeNoise = changeNoise;
condition.noChangeNoise = noChangeNoise;
condition.changeComponentCoherence = changeSignal + changeNoise;
condition.noChangeComponentCoherence = noChangeSignal + noChangeNoise;
condition.changeSignalCoherence = changeSignal;
condition.noChangeSignalCoherence = noChangeSignal;
condition.preStepComponentCoherence = repmat(pedestalByTrial,1,3);
condition.experimentalCorrect = logical(T.correct);
condition.nExcludedNonfinite = nExcluded;
condition.pedestalPC = idqPedestalPC(:)';
condition.definition = [ ...
  "IDQ drift-relative rectangular step predictors are signed deviations " + ...
  "about a positive coherence pedestal in all three directions on both " + ...
  "patches; stepCoh is added to the changed-patch drift pedestal"];
end

function conditions = makeIDR(S)
if ~isfield(S,'sideGainData')
  error('makeIDRIDQMTTrialInputs:MissingIDRSideGainData', ...
    'IDR inventory has no sideGainData. Rerun analyzeIDRSideGains.');
end
D = S.sideGainData; T = D.trialTable;
vars = {'sessionIndex','trialIdx','probeDirDeg','stepCohPC', ...
  'stepDeltaCohPC','nYokedProbeStreams','correct'};
requireTableVariables(T,vars,'IDR');
required = {'stepFrames','changePrefNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial','changeProbeCandidateNoiseByFrameTrial', ...
  'noChangeProbeCandidateNoiseByFrameTrial'};
if ~all(isfield(D,required))
  error('makeIDRIDQMTTrialInputs:MissingIDRNoise', ...
    'IDR sideGainData lacks required per-candidate noise matrices.');
end

prefC = stepMean(D.changePrefNoiseByFrameTrial,D.stepFrames);
prefN = stepMean(D.noChangePrefNoiseByFrameTrial,D.stepFrames);
probeC = stepMean(D.changeProbeCandidateNoiseByFrameTrial,D.stepFrames);
probeN = stepMean(D.noChangeProbeCandidateNoiseByFrameTrial,D.stepFrames);
offsets = unique(double(T.probeDirDeg),'sorted');
offsets = offsets(isfinite(offsets) & offsets ~= 1);
conditions = repmat(baseCondition("IDR",nan,[],strings(0,1),[],[]), ...
  numel(offsets),1);
for k = 1:numel(offsets)
  theta = offsets(k); use = double(T.probeDirDeg)==theta;
  n = sum(use); nStreams = unique(double(T.nYokedProbeStreams(use)));
  if numel(nStreams)~=1 || ~ismember(nStreams,[1 2])
    error('makeIDRIDQMTTrialInputs:BadYoking', ...
      'Offset %.6g has inconsistent probe-stream counts.',theta);
  end
  if nStreams==2
    directions = [0 theta -theta];
    roles = ["primary","probePlus","probeMinus"];
    cNoise = [prefC(use),probeC(use),probeC(use)];
    nNoise = [prefN(use),probeN(use),probeN(use)];
  else
    directions = [0 180]; roles = ["primary","probe"];
    cNoise = [prefC(use),probeC(use)];
    nNoise = [prefN(use),probeN(use)];
  end
  post = double(T.stepCohPC(use));
  base = post-double(T.stepDeltaCohPC(use));
  cSignal = [post zeros(n,numel(directions)-1)];
  nSignal = [base zeros(n,numel(directions)-1)];
  finiteRows=all(isfinite(cNoise),2)&all(isfinite(nNoise),2)& ...
    all(isfinite(cSignal),2)&all(isfinite(nSignal),2);
  nExcluded=sum(~finiteRows);
  cNoise=cNoise(finiteRows,:);nNoise=nNoise(finiteRows,:);
  cSignal=cSignal(finiteRows,:);nSignal=nSignal(finiteRows,:);
  C = baseCondition("IDR",theta,directions,roles, ...
    T.sessionIndex(use,:),T.trialIdx(use,:));
  C.sessionIndex=C.sessionIndex(finiteRows);
  C.trialIndex=C.trialIndex(finiteRows);
  C.changeNoise = cNoise; C.noChangeNoise = nNoise;
  C.changeComponentCoherence = cSignal+cNoise;
  C.noChangeComponentCoherence = nSignal+nNoise;
  C.changeSignalCoherence = cSignal;
  C.noChangeSignalCoherence = nSignal;
  C.preStepComponentCoherence = nSignal;
  correct=logical(T.correct(use));
  C.experimentalCorrect=correct(finiteRows);
  C.nExcludedNonfinite = nExcluded;
  C.definition = [ ...
    "IDR drift-relative rectangular step predictors; probeCandidate is " + ...
    "replicated at +theta and -theta only for two-stream yoked probes; " + ...
    "changed and unchanged preferred signals use stepCohPC and " + ...
    "stepCohPC-stepDeltaCohPC"];
  conditions(k) = C;
end
end

function C = baseCondition(dataset,offset,directions,roles,sessionIndex,trialIndex)
C = struct('dataset',dataset,'offsetDeg',offset, ...
  'componentDirectionsDeg',double(directions(:)'), ...
  'componentRoles',string(roles(:)'), ...
  'sessionIndex',double(sessionIndex(:)), ...
  'trialIndex',double(trialIndex(:)), ...
  'changeNoise',[],'noChangeNoise',[], ...
  'changeComponentCoherence',[],'noChangeComponentCoherence',[], ...
  'changeSignalCoherence',[],'noChangeSignalCoherence',[], ...
  'preStepComponentCoherence',[], ...
  'experimentalCorrect',false(0,1), ...
  'nExcludedNonfinite',0,'pedestalPC',[],'definition',"");
end

function pedestal=resolveIDQPedestal(values,sessionIndex)
values=double(values(:)); sessionIndex=double(sessionIndex(:));
if isscalar(values)
  pedestal=repmat(values,numel(sessionIndex),1);
elseif all(sessionIndex==fix(sessionIndex)) && all(sessionIndex>=1) && ...
    numel(values)>=max(sessionIndex)
  pedestal=values(sessionIndex);
else
  error('makeIDRIDQMTTrialInputs:BadIDQPedestalSize', ...
    ['IDQPedestalPC must be scalar or indexed by every positive integer ' ...
     'sessionIndex in the IDQ trial table.']);
end
end

function x = stepMean(x,frames)
x = mean(double(x(frames,:)),1,'omitnan')';
end

function X = tableColumns(T,names)
X = nan(height(T),numel(names));
for k=1:numel(names), X(:,k)=double(T.(names{k})); end
end

function requireTableVariables(T,names,label)
missing = setdiff(names,T.Properties.VariableNames);
if ~isempty(missing)
  error('makeIDRIDQMTTrialInputs:MissingTableVariables', ...
    '%s trial table is missing: %s',label,strjoin(missing,', '));
end
end
