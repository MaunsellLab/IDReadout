function [summary, figures] = plotIDLAcrossSessionBehavior(varargin)
% plotIDLAcrossSessionBehavior  Descriptive IDL functions and noise kernels.
% TwoAllocationMode=false selects the original style; true selects 100:0/50:50.
%
% [summary, figures] = plotIDLAcrossSessionBehavior()
% [summary, figures] = plotIDLAcrossSessionBehavior('Animal','Meetz')
%
% Reads IDL/Data/SessionAnalysis using selectAnalysisFiles. Produces two
% docked figure tabs per group (Meetz, Neesha, and, when both are selected,
% Combined): raw and aligned psychometric functions; empirical kernels for
% seven signal fractions. Also exports PDFs to IDL/Plots/AcrossSessions.
% Set 'ExportPDF',false to skip the files, or 'Docked',false for independent
% figure windows. Other name-value pairs are passed to selectAnalysisFiles.
%
% Alignment uses one multiplicative step shift per session, estimated as
% the geometric mean of its four fitted 75%-correct thresholds. It cannot
% create or erase differences among allocations within a session. Because
% these choices also determine the alignment, aligned curves are displays,
% not independent validation of slope or threshold differences. Displayed
% group functions fit a separate Weibull threshold and slope per allocation;
% they use a common group lapse from fitIDLSessionPsychometrics.
%
% The gain regression uses one predictor per stream and patch: the realized
% coherence at all 30 video frames, integrated with the session's LOO model
% kernel. Descriptive kernels are correct-minus-error mean coherence at
% every video frame. Both estimates remove session-by-step noise means.
% Scalar linear-probability coefficients are retained as integratedGain.
% Displayed relative gains divide by the observed 100%-Change coefficient
% in the same animal/group; they are not psychometric-fit gain parameters.
% Prestep and step noise are both reconstructed at every video frame.
% The allocation panel also projects per-trial model probabilities from
% Data/AcrossSessions/ModelRecovery.mat through the scalar regression.
% Run fitIDLModelRecoveryAcrossSessions on the current sessions first, or
% pass its returned result with 'ModelResult',result. No synthetic choices
% are needed for this conditional expected regression coefficient.
% Four prediction curves are displayed: L1, L2, Max, and Square (sum of
% squared integrated stream responses, without the L2 square root).
% Re-run fitIDLModelRecoveryAcrossSessions after installing this update;
% older ModelRecovery files lack the required pSquare predictions.
% Noise-free trials remain in psychometric displays, but trials without
% recorded noise anywhere in the trial do not enter noise kernels/gains.
% Unidentifiable gain coefficients are NaN, not zero. If there is no
% estimable observed 100%-Change reference, relative gains are unavailable.
% Primary analysis excludes sessions without recorded step noise, matching
% fitIDLModelRecoveryAcrossSessions. IncludeNoiseFreeSessions=true restores
% those sessions for secondary behavioral displays (not noise gains).

p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;
p.addParameter('SessionFolder','',@(x) ischar(x)||isstring(x));
p.addParameter('TwoAllocationMode',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('ExportPDF', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Docked', true, @(x) islogical(x) && isscalar(x));
p.addParameter('Visible', 'on', @(x) any(strcmpi(string(x),["on" "off"])));
p.addParameter('ModelResult', [], @(x) isempty(x) || isstruct(x));
p.addParameter('ModelRecoveryPath', '', @(x) ischar(x) || isstring(x));
p.addParameter('IncludeNoiseFreeSessions', false, ...
  @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;
tag = idlAllocationTag(opts.TwoAllocationMode);

root = domainFolder(mfilename('fullpath'));
dataFolder=char(string(opts.SessionFolder));
if isempty(dataFolder), dataFolder = fullfile(root,'Data','SessionAnalysis'); end
assert(isfolder(dataFolder), 'plotIDLAcrossSessionBehavior:MissingData', ...
  'Cannot find the IDL SessionAnalysis folder: %s', dataFolder);
selectionArguments = unmatchedPairs(p.Unmatched);
sourceList = selectAnalysisFiles(dataFolder,'TwoAllocationMode',opts.TwoAllocationMode,selectionArguments{:});
if ischar(sourceList), sourceList = {sourceList}; end
if isstring(sourceList), sourceList = cellstr(sourceList(:)); end
sourceList = sourceList(:);
assert(~isempty(sourceList), 'plotIDLAcrossSessionBehavior:NoSessions', ...
  'No sessionAnalysis files were selected.');

sessions = cell(numel(sourceList),1);
sessionTables = cell(numel(sourceList),1);
allocationPC = [100 75 62.5 50];
if opts.TwoAllocationMode, allocationPC=[100 50]; end
signalFractionPC = sort(unique([allocationPC 100-allocationPC]));
animal = strings(numel(sourceList),1);
sessionID = strings(numel(sourceList),1);
nTrials = zeros(numel(sourceList),1);
thresholds = nan(numel(sourceList),4);
beta = nan(numel(sourceList),1);
lapse = nan(numel(sourceList),1);
scale = nan(numel(sourceList),1);
selected = true(numel(sourceList),1);
for s = 1:numel(sourceList)
  loaded = load(sourceList{s}, 'sessionAnalysis');
  assert(isfield(loaded,'sessionAnalysis'), ...
    'plotIDLAcrossSessionBehavior:MissingAnalysis', ...
    '%s lacks sessionAnalysis.', sourceList{s});
  A = loaded.sessionAnalysis;
  idlAssertAnalysisMode(A,opts.TwoAllocationMode);
  hasNoise = any(abs(double(A.changeNoiseByFrameTrial))>1e-8,'all') || ...
    any(abs(double(A.noChangeNoiseByFrameTrial))>1e-8,'all');
  if ~opts.IncludeNoiseFreeSessions && ~hasNoise
    selected(s) = false;
    continue
  end
  A = ensureFullNoise(A,sourceList{s});
  validateSession(A, sourceList{s});
  sessions{s} = A;
  animal(s) = string(A.sessionHeader.animal);
  sessionID(s) = string(A.sourceSessionID);
  T = A.trialTable;
  nTrials(s) = height(T);
  fit = fitIDLSessionPsychometrics(T);
  assert(all(isfinite(fit.thresholdTable.threshold)) && ...
    all(fit.thresholdTable.threshold > 0), ...
    'plotIDLAcrossSessionBehavior:ThresholdFit', ...
    'Cannot align session %s: nonpositive or nonfinite threshold.', sessionID(s));
  [~,columns] = ismember(allocationPC,[100 75 62.5 50]);
  thresholds(s,columns) = fit.thresholdTable.threshold(:)';
  beta(s) = fit.betaWeibull;
  lapse(s) = fit.lapse;
  scale(s) = exp(mean(log(thresholds(s,:)),'omitnan'));
  T.animal = repmat(animal(s),height(T),1);
  T.sourceSessionID = repmat(sessionID(s),height(T),1);
  T.alignedStep = double(T.totalStepCohPC) ./ scale(s);
  sessionTables{s} = T;
end

excludedNoiseFreeSourceFiles = string(sourceList(~selected));
sourceList = sourceList(selected);
sessions = sessions(selected);
sessionTables = sessionTables(selected);
animal = animal(selected); sessionID = sessionID(selected);
nTrials = nTrials(selected); thresholds = thresholds(selected,:);
beta = beta(selected); lapse = lapse(selected); scale = scale(selected);
assert(~isempty(sourceList), 'plotIDLAcrossSessionBehavior:NoNoisySessions', ...
  ['No sessions with recorded step noise were selected. Set ' ...
   'IncludeNoiseFreeSessions=true for secondary noise-free displays.']);

modelResult = opts.ModelResult;
if isempty(modelResult)
  modelPath = string(opts.ModelRecoveryPath);
  if strlength(modelPath)==0
    modelPath = fullfile(root,'Data','AcrossSessions',['ModelRecovery' tag '.mat']);
  end
  assert(isfile(modelPath), 'plotIDLAcrossSessionBehavior:MissingModelRecovery', ...
    ['Cannot find %s. Run fitIDLModelRecoveryAcrossSessions(''Animal'',''All'') ' ...
    'on the current sessions first, or pass ''ModelResult'',result.'],modelPath);
  loaded = load(modelPath,'result');
  assert(isfield(loaded,'result'), ...
    'plotIDLAcrossSessionBehavior:MissingModelResult', ...
    '%s does not contain a result struct.',modelPath);
  modelResult = loaded.result;
end
assert(isfield(modelResult,'twoAllocationMode') && ...
  modelResult.twoAllocationMode==opts.TwoAllocationMode, ...
  'plotIDLAcrossSessionBehavior:ModelStyle','Model result has a different allocation style.');
modelProbabilities = matchModelProbabilities(modelResult,sessions);
sessionWeights = matchModelWeights(modelResult,sessions);

trials = vertcat(sessionTables{:});
animals = unique(animal,'stable');
groups = [animals; "Combined"];
if numel(animals) == 1, groups = animals; end
groupData = struct([]);
for g = 1:numel(groups)
  members = animal == groups(g);
  if groups(g) == "Combined", members(:) = true; end
  Tg = vertcat(sessionTables{members});
  groupData(g).name = groups(g);
  groupData(g).twoAllocationMode=opts.TwoAllocationMode;
  groupData(g).sessionIndex = find(members);
  groupData(g).trials = Tg;
  groupData(g).rawFit = fitGroupFunctions(Tg,'totalStepCohPC',allocationPC);
  groupData(g).alignedFit = fitGroupFunctions(Tg,'alignedStep',allocationPC);
  groupData(g).kernel = estimateKernels( ...
    sessions,members,signalFractionPC,modelProbabilities,sessionWeights);
end

sessionThresholdTable = table(sessionID,animal,nTrials, ...
  thresholds(:,1),thresholds(:,2),thresholds(:,3),thresholds(:,4), ...
  beta,lapse,scale, ...
  'VariableNames',{'sourceSessionID','animal','nTrials', ...
  'threshold100','threshold75','threshold62_5','threshold50', ...
  'betaWeibull','lapse','geometricMeanThreshold'});

rawMax = max(double(trials.totalStepCohPC));
alignedMax = max(double(trials.alignedStep));
rawLimit = [0 max(5,ceil(rawMax/5)*5)];
alignedLimit = [0 max(1,ceil(alignedMax*2)/2)];
kernelScale = max(arrayfun(@(g) max(abs(g.kernel.coefficients),[], ...
  'all','omitnan'),groupData));
kernelScale = max([kernelScale eps],[],'omitnan');
summaryScale = max(arrayfun(@(g) max(abs([ ...
  g.kernel.relativeGain(:); ...
  g.kernel.predictedRelativeGain(:)]),[],'omitnan'),groupData));
summaryScale = max([summaryScale eps],[],'omitnan');

figures = struct();
for g = 1:numel(groups)
  key = matlab.lang.makeValidName(char(groups(g)));
  figF = newFigure("IDL psychometric functions | " + groups(g) + string(tag),opts);
  plotFunctions(figF,groupData(g),groupData,sessionThresholdTable, ...
    allocationPC,rawLimit,alignedLimit);
  figK = newFigure("IDL behavioral kernels | " + groups(g) + string(tag),opts);
  plotKernels(figK,groupData(g),signalFractionPC,kernelScale,summaryScale);
  figures.(key) = struct('functions',figF,'kernels',figK);
end

if opts.ExportPDF
  plotFolder = validFolder(fullfile(root,'Plots','AcrossSessions'));
  for g = 1:numel(groups)
    key = matlab.lang.makeValidName(char(groups(g)));
    exportgraphics(figures.(key).functions, ...
      fullfile(plotFolder,['IDLFunctions_' key tag '.pdf']),'ContentType','vector');
    exportgraphics(figures.(key).kernels, ...
      fullfile(plotFolder,['IDLKernels_' key tag '.pdf']),'ContentType','vector');
  end
end

summary = struct();
summary.twoAllocationMode = opts.TwoAllocationMode;
summary.sourceFiles = string(sourceList);
summary.includeNoiseFreeSessions = opts.IncludeNoiseFreeSessions;
summary.excludedNoiseFreeSourceFiles = excludedNoiseFreeSourceFiles;
summary.sessionThresholdTable = sessionThresholdTable;
summary.groups = rmfield(groupData,'trials');
summary.signalFractionPC = signalFractionPC;
summary.kernelTimeCentersMS = groupData(1).kernel.timeCentersMS;
summary.kernelUnits = 'mean coherence noise on correct minus error trials (percentage points)';
summary.gainUnits = 'integratedGain: change in P(correct) per coherence percentage point; relativeGain: integratedGain / observed 100%-Change integratedGain';
summary.alignment = 'total step divided by session geometric mean of enabled-allocation 75% thresholds';
summary.modelPredictionSource = 'ModelRecovery per-trial probabilities; fixed fitted kernels and sigmas';
if isfield(modelResult,'sigmaFitMode')
  summary.sigmaFitMode = modelResult.sigmaFitMode;
  summary.modelScoreDescription = modelResult.scoreDescription;
end
end

%% ------------------------------------------------------------------------
function pairs = unmatchedPairs(unmatched)
names = fieldnames(unmatched);
pairs = cell(1,2*numel(names));
for k = 1:numel(names)
  pairs{2*k-1} = names{k};
  pairs{2*k} = unmatched.(names{k});
end
end

%% ------------------------------------------------------------------------
function probabilities = matchModelProbabilities(modelResult,sessions)
assert(isscalar(modelResult) && isfield(modelResult,'predictionTable') && ...
  isfield(modelResult,'modelNames') && ...
  isequal(string(modelResult.modelNames(:))',["L1" "L2" "Linf" "Square"]), ...
  'plotIDLAcrossSessionBehavior:ModelResultFormat', ...
  ['ModelResult must contain L1, L2, Linf, and Square. Re-run the updated ' ...
   'fitIDLModelRecoveryAcrossSessions on the current sessions.']);
P = modelResult.predictionTable;
required = {'sourceSessionID','sourceTrialNumber','correct', ...
  'unsignedAllocPC','totalStepCohPC','pL1','pL2','pLinf','pSquare'};
assert(istable(P) && all(ismember(required,P.Properties.VariableNames)), ...
  'plotIDLAcrossSessionBehavior:PredictionFields', ...
  ['The saved predictionTable is missing required trial fields, including ' ...
   'pSquare. Re-run the updated fitIDLModelRecoveryAcrossSessions.']);
probabilities = cell(numel(sessions),1);
for s = 1:numel(sessions)
  A = sessions{s};
  T = A.trialTable;
  rows = string(P.sourceSessionID)==string(A.sourceSessionID);
  Q = P(rows,:);
  [found,order] = ismember(double(T.sourceTrialNumber(:)), ...
    double(Q.sourceTrialNumber(:)));
  assert(height(Q)==height(T) && all(found) && ...
    numel(unique(Q.sourceTrialNumber))==height(T), ...
    'plotIDLAcrossSessionBehavior:ModelTrialMismatch', ...
    ['ModelRecovery lacks an exact trial match for %s. Re-run ' ...
    'fitIDLModelRecoveryAcrossSessions on the current sessions.'], ...
    string(A.sourceSessionID));
  Q = Q(order,:);
  assert(all(double(Q.correct)==double(T.correct)) && ...
    all(abs(double(Q.unsignedAllocPC)-double(T.unsignedAllocPC))<1e-5) && ...
    all(abs(double(Q.totalStepCohPC)-double(T.totalStepCohPC))<1e-5), ...
    'plotIDLAcrossSessionBehavior:StaleModelRecovery', ...
    'Choices or steps changed for %s; re-run fitIDLModelRecoveryAcrossSessions.', ...
    string(A.sourceSessionID));
  if ismember('prestepCohPC',Q.Properties.VariableNames)
    assert(all(abs(double(Q.prestepCohPC)-double(T.prestepCohPC))<1e-5), ...
      'plotIDLAcrossSessionBehavior:StalePedestal', ...
      'Pedestals changed for %s; re-run fitIDLModelRecoveryAcrossSessions.', ...
      string(A.sourceSessionID));
  end
  if ismember('noiseCohPC',Q.Properties.VariableNames) && ...
      ismember('noiseCohPC',T.Properties.VariableNames)
    assert(all(abs(double(Q.noiseCohPC)-double(T.noiseCohPC))<1e-5), ...
      'plotIDLAcrossSessionBehavior:StaleNoise', ...
      'Noise settings changed for %s; re-run fitIDLModelRecoveryAcrossSessions.', ...
      string(A.sourceSessionID));
  end
  if ismember('hasRecordedStepNoise',Q.Properties.VariableNames)
    C = any(any(abs(double(A.changeNoiseByFrameTrial))>1e-8,1),3);
    N = any(any(abs(double(A.noChangeNoiseByFrameTrial))>1e-8,1),3);
    assert(isequal(logical(Q.hasRecordedStepNoise),reshape(C | N,[],1)), ...
      'plotIDLAcrossSessionBehavior:StaleRecordedNoise', ...
      'Recorded noise changed for %s; re-run fitIDLModelRecoveryAcrossSessions.', ...
      string(A.sourceSessionID));
  end
  probabilities{s} = double(Q{:,{'pL1','pL2','pLinf','pSquare'}});
  assert(all(isfinite(probabilities{s}),'all') && ...
    all(probabilities{s}>=0,'all') && all(probabilities{s}<=1,'all'), ...
    'plotIDLAcrossSessionBehavior:InvalidProbabilities', ...
    'Nonfinite or out-of-range model probabilities in %s.', ...
    string(A.sourceSessionID));
end
end

%% ------------------------------------------------------------------------
function weights = matchModelWeights(modelResult,sessions)
assert(isfield(modelResult,'sessionTable') && ...
  ismember('sourceSessionID',modelResult.sessionTable.Properties.VariableNames) && ...
  isfield(modelResult,'kernelByFrameSession'), ...
  'plotIDLAcrossSessionBehavior:MissingModelKernels', ...
  'ModelResult must contain sessionTable and kernelByFrameSession.');
ids = string(modelResult.sessionTable.sourceSessionID);
W = double(modelResult.kernelByFrameSession);
assert(size(W,2)==numel(ids), ...
  'plotIDLAcrossSessionBehavior:KernelSessionCount', ...
  'Kernel and session counts differ in ModelResult.');
weights = cell(numel(sessions),1);
for s = 1:numel(sessions)
  index = find(ids==string(sessions{s}.sourceSessionID));
  assert(numel(index)==1 && size(W,1)==numel(sessions{s}.stepTimeMS), ...
    'plotIDLAcrossSessionBehavior:KernelMismatch', ...
    'Missing or incompatible temporal kernel for %s.', ...
    string(sessions{s}.sourceSessionID));
  weights{s} = W(:,index);
  assert(all(isfinite(weights{s})) && abs(sum(abs(weights{s}))-1)<1e-5, ...
    'plotIDLAcrossSessionBehavior:KernelNormalization', ...
    'Expected a finite, unit-L1 temporal kernel for %s.', ...
    string(sessions{s}.sourceSessionID));
end
end

%% ------------------------------------------------------------------------
function validateSession(A,source)
required = {'trialTable','sourceSessionID','sessionHeader','stepTimeMS', ...
  'changeNoiseByFrameTrial','noChangeNoiseByFrameTrial'};
assert(all(isfield(A,required)), ...
  'plotIDLAcrossSessionBehavior:MissingFields', ...
  'Missing sessionAnalysis fields in %s.',source);
T = A.trialTable;
assert(all(ismember({'correct','totalStepCohPC','unsignedAllocPC','alloc0'}, ...
  T.Properties.VariableNames)), ...
  'plotIDLAcrossSessionBehavior:MissingTrialFields', ...
  'Missing trialTable fields in %s.',source);
n = height(T);
nf = numel(A.stepTimeMS);
assert(isequal(size(A.changeNoiseByFrameTrial),[nf n 2]) && ...
  isequal(size(A.noChangeNoiseByFrameTrial),[nf n 2]), ...
  'plotIDLAcrossSessionBehavior:NoiseShape', ...
  'Noise arrays do not match trial counts in %s.',source);
expectedFrames=round(double(A.sessionHeader.stepMS)*double(A.sessionHeader.frameRateHz)/1000);
assert(nf==expectedFrames, 'plotIDLAcrossSessionBehavior:TimeGrid', ...
  'Expected %d video frames in %s; found %d.',expectedFrames,source,nf);
assert(isfield(A,'fullTimeMS') && ...
  isequal(size(A.fullChangeNoiseByFrameTrial),[numel(A.fullTimeMS) n 2]) && ...
  isequal(size(A.fullNoChangeNoiseByFrameTrial),[numel(A.fullTimeMS) n 2]), ...
  'plotIDLAcrossSessionBehavior:FullNoiseShape', ...
  'Full-trial noise arrays do not match trials in %s.',source);
idlAssertAnalysisMode(A,logical(A.sessionHeader.twoAllocationMode));
end

%% ------------------------------------------------------------------------
function A = ensureFullNoise(A,compactPath)
if all(isfield(A,{'fullTimeMS','fullChangeNoiseByFrameTrial', ...
    'fullNoChangeNoiseByFrameTrial'}))
  return
end
assert(isfield(A,'sourceSessionPath') && ...
  isfile(A.sourceSessionPath) && ...
  ismember('sourceTrialIndex',A.trialTable.Properties.VariableNames), ...
  'plotIDLAcrossSessionBehavior:MissingFullSession', ...
  ['Full prestep noise is absent from %s. Its sourceSessionPath must ' ...
   'point to the original FullSession file, or regenerate this compact ' ...
   'session with the updated makeIDLSessionAnalysis.'],compactPath);
source = load(A.sourceSessionPath,'trials');
frameMS = 1000/double(A.sessionHeader.frameRateHz);
nFrames = round((double(A.sessionHeader.preStepMS)+ ...
  double(A.sessionHeader.stepMS))/frameMS);
times = (0:nFrames-1)'*frameMS;
n = height(A.trialTable);
A.fullTimeMS = times-double(A.sessionHeader.preStepMS);
A.fullChangeNoiseByFrameTrial = zeros(nFrames,n,2,'single');
A.fullNoChangeNoiseByFrameTrial = zeros(nFrames,n,2,'single');
for k = 1:n
  tr = source.trials{double(A.trialTable.sourceTrialIndex(k))};
  A.fullChangeNoiseByFrameTrial(:,k,:) = frameNoise(tr,'change',times);
  A.fullNoChangeNoiseByFrameTrial(:,k,:) = frameNoise(tr,'noChange',times);
end
step = A.fullTimeMS>=0;
assert(isequal(A.changeNoiseByFrameTrial, ...
  A.fullChangeNoiseByFrameTrial(step,:,:)) && ...
  isequal(A.noChangeNoiseByFrameTrial, ...
  A.fullNoChangeNoiseByFrameTrial(step,:,:)), ...
  'plotIDLAcrossSessionBehavior:NoiseReconstruction', ...
  'Reconstructed step noise differs from the compact analysis in %s.', ...
  compactPath);
end

%% ------------------------------------------------------------------------
function values = frameNoise(tr,prefix,times)
eventTime = double(tr.([prefix 'TimesMS']).data(:));
values = zeros(numel(times),1,2,'single');
for stream = 1:2
  field = sprintf('%sNoise%dCohsPC',prefix,stream-1);
  coherence = double(tr.(field).data(:));
  assert(numel(coherence)==numel(eventTime), ...
    'plotIDLAcrossSessionBehavior:NoiseEvents', ...
    '%s has mismatched noise event times and coherence.',prefix);
  for f = 1:numel(times)
    index = find(eventTime<=times(f)+1e-6,1,'last');
    assert(~isempty(index), ...
      'plotIDLAcrossSessionBehavior:MissingNoiseEvent', ...
      'Missing %s coherence at frame %d.',prefix,f);
    values(f,1,stream) = single(coherence(index));
  end
end
end

%% ------------------------------------------------------------------------
function F = fitGroupFunctions(T,stepField,allocationPC)
% One common lapse stabilizes the per-allocation descriptive slope fits.
Tfit = T;
Tfit.totalStepCohPC = double(T.(stepField));
base = fitIDLSessionPsychometrics(Tfit);
F = struct('allocationPC',allocationPC(:), ...
  'threshold',nan(numel(allocationPC),1),'beta',nan(numel(allocationPC),1), ...
  'alpha',nan(numel(allocationPC),1),'lapse',base.lapse,'nTrials',nan(numel(allocationPC),1));
for k = 1:numel(allocationPC)
  sel = abs(double(T.unsignedAllocPC)-allocationPC(k))<1e-6;
  stepValues = double(T.(stepField));
  x = stepValues(sel);
  y = double(T.correct(sel));
  F.nTrials(k) = numel(x);
  start = log([base.thresholdTable.threshold(k);base.betaWeibull]);
  objective = @(z) curveNLL(z,x,y,base.lapse);
  options = optimset('Display','off','MaxIter',3000,'MaxFunEvals',3000);
  z = fminsearch(objective,start,options);
  F.threshold(k) = exp(z(1));
  F.beta(k) = exp(z(2));
  F.alpha(k) = idqWeibullAlphaForThreshold(0.75, ...
    F.threshold(k),F.beta(k),F.lapse);
end
end

%% ------------------------------------------------------------------------
function nll = curveNLL(z,x,y,lapse)
theta = exp(z(1));
beta = exp(z(2));
if ~isfinite(theta) || ~isfinite(beta) || theta < 1e-5 || beta < 0.05 || beta > 30
  nll = 1e12;
  return
end
alpha = idqWeibullAlphaForThreshold(0.75,theta,beta,lapse);
prob = idqWeibullP(x,alpha,beta,lapse);
prob = min(max(prob,1e-12),1-1e-12);
nll = -sum(y.*log(prob)+(1-y).*log1p(-prob));
end

%% ------------------------------------------------------------------------
function K = estimateKernels(sessions,members,signalFractionPC,modelProbabilities,sessionWeights)
% At each allocation, fit stream/patch integrated predictors jointly.
% Framewise traces are conditional correct-minus-error coherence means.
unsigned = sort(unique(sessions{find(members,1)}.trialTable.unsignedAllocPC),'descend')';
nFrames = numel(sessions{1}.fullTimeMS);
stepFrames = sessions{1}.fullTimeMS >= 0;
nStepFrames = sum(stepFrames);
nResponses = 1+size(modelProbabilities{find(members,1)},2);
frameSet = nan(numel(signalFractionPC),nFrames,2,nResponses);
gainSet = nan(numel(signalFractionPC),2,nResponses);
counts = zeros(numel(signalFractionPC),1);
for a = 1:numel(unsigned)
  highFraction = unsigned(a);
  equal = highFraction == 50;
  if equal, nStreams = 1; else, nStreams = 2; end
  nPredictors = 2*nStreams;
  sumCorrect = zeros(nPredictors,nResponses,nFrames);
  sumError = zeros(nPredictors,nResponses,nFrames);
  countCorrect = zeros(1,nResponses);
  countError = zeros(1,nResponses);
  XXgain = zeros(nPredictors);
  XyGain = zeros(nPredictors,nResponses);
  nTrials = 0;
  for s = find(members(:))'
    A = sessions{s}; T = A.trialTable;
    rows = find(abs(double(T.unsignedAllocPC)-highFraction)<1e-6);
    % A genuinely noise-free trial has no information about noise gain.
    % Keep it in psychometrics and likelihood fits, not in noise statistics.
    C = any(any(abs(double(A.fullChangeNoiseByFrameTrial(:,rows,:)))>1e-8,1),3);
    N = any(any(abs(double(A.fullNoChangeNoiseByFrameTrial(:,rows,:)))>1e-8,1),3);
    rows = rows(reshape(C | N,[],1));
    if isempty(rows), continue; end
    n = numel(rows); nTrials = nTrials+n;
    xFrame = zeros(n,nPredictors,nFrames);
    highIs0 = double(T.alloc0(rows))>double(T.alloc1(rows));
    for side = 1:2
      if side == 1, noise = A.fullChangeNoiseByFrameTrial;
      else, noise = A.fullNoChangeNoiseByFrameTrial;
      end
      column = (side-1)*nStreams;
      for f = 1:nFrames
        n0 = reshape(double(noise(f,rows,1)),[],1);
        n1 = reshape(double(noise(f,rows,2)),[],1);
        if equal
          xFrame(:,column+1,f) = n0+n1;
        else
          high = n0; low = n1;
          high(~highIs0) = n1(~highIs0);
          low(~highIs0) = n0(~highIs0);
          xFrame(:,column+1,f) = high;
          xFrame(:,column+2,f) = low;
        end
      end
    end
    probability = modelProbabilities{s};
    yRaw = [double(T.correct(rows)),probability(rows,:)];
    y = yRaw;
    step = double(T.totalStepCohPC(rows));
    assert(all(isfinite(xFrame),'all') && all(isfinite(yRaw),'all') && ...
      all(isfinite(step)), ...
      'plotIDLAcrossSessionBehavior:NonfiniteTrial', ...
      'Nonfinite choice, step, or noise in %s.',string(A.sourceSessionID));
    % Remove every session-by-step intercept from choices and predictors.
    for level = unique(step(:))'
      within = step == level;
      y(within,:) = y(within,:)-mean(y(within,:),1);
      xFrame(within,:,:) = xFrame(within,:,:)-mean(xFrame(within,:,:),1);
    end
    w = sessionWeights{s};
    assert(isequal(A.fullTimeMS>=0,stepFrames), ...
      'plotIDLAcrossSessionBehavior:FullTimeGrid', ...
      'Inconsistent full trial time grid across sessions.');
    xStep = xFrame(:,:,stepFrames);
    xGain = reshape(reshape(xStep,n*nPredictors,nStepFrames)*w,n,nPredictors);
    XXgain = XXgain+xGain'*xGain;
    XyGain = XyGain+xGain'*y;
    for f = 1:nFrames
      x = xFrame(:,:,f);
      sumCorrect(:,:,f) = sumCorrect(:,:,f)+x'*yRaw;
      sumError(:,:,f) = sumError(:,:,f)+x'*(1-yRaw);
    end
    countCorrect = countCorrect+sum(yRaw,1);
    countError = countError+sum(1-yRaw,1);
  end
  fittedGain = nan(nPredictors,nResponses);
  active = diag(XXgain)>1e-10;
  if any(active)
    assert(nTrials>sum(active) && rcond(XXgain(active,active))>1e-10, ...
      'plotIDLAcrossSessionBehavior:GainRank', ...
      'Too few independent integrated noise samples at %g:%g.', ...
      highFraction,100-highFraction);
    fittedGain(active,:) = XXgain(active,active)\XyGain(active,:);
  end
  fittedFrame = sumCorrect./reshape(countCorrect,1,nResponses,1) - ...
    sumError./reshape(countError,1,nResponses,1);
  if equal, fittedFrame = fittedFrame/2; end % two exchangeable streams
  if equal
    index = find(signalFractionPC==50);
    for side = 1:2
      col = side;
      gainSet(index,side,:) = reshape(fittedGain(col,:),1,1,nResponses);
      frameSet(index,:,side,:) = permute(fittedFrame(col,:,:),[4 3 1 2]);
    end
    counts(index) = nTrials;
  else
    highIndex = find(abs(signalFractionPC-highFraction)<1e-6);
    lowIndex = find(abs(signalFractionPC-(100-highFraction))<1e-6);
    for side = 1:2
      highCol = (side-1)*2+1;
      lowCol = highCol+1;
      gainSet(highIndex,side,:) = reshape(fittedGain(highCol,:),1,1,nResponses);
      gainSet(lowIndex,side,:) = reshape(fittedGain(lowCol,:),1,1,nResponses);
      frameSet(highIndex,:,side,:) = permute(fittedFrame(highCol,:,:),[4 3 1 2]);
      frameSet(lowIndex,:,side,:) = permute(fittedFrame(lowCol,:,:),[4 3 1 2]);
    end
    counts([highIndex lowIndex]) = nTrials;
  end
end
reference = gainSet(signalFractionPC==100,1,1);
referenceAvailable = isfinite(reference) && reference>1e-8;
relativeGain = nan(size(gainSet(:,:,1)));
predictedRelativeGain = nan(size(gainSet(:,:,2:nResponses)));
if referenceAvailable
  relativeGain = gainSet(:,:,1)/reference;
  predictedRelativeGain = gainSet(:,:,2:nResponses)/reference;
end
K = struct('signalFractionPC',signalFractionPC(:), ...
  'timeCentersMS',double(sessions{1}.fullTimeMS(:))', ...
  'preStepMS',double(sessions{1}.sessionHeader.preStepMS), ...
  'coefficients',frameSet(:,:,:,1), ...
  'predictedCoefficients',frameSet(:,:,:,2:nResponses), ...
  'integratedGain',gainSet(:,:,1), ...
  'predictedIntegratedGain',gainSet(:,:,2:nResponses), ...
  'referenceGain',reference, ...
  'referenceGainAvailable',referenceAvailable, ...
  'relativeGain',relativeGain, ...
  'predictedRelativeGain',predictedRelativeGain, ...
  'nTrials',counts, 'patchSource',{{'Change','No change'}}, ...
  'predictedModelNames',["L1" "L2" "Max" "Square"]);
K.timeAveragedCoefficients = squeeze(mean(K.coefficients,2));
K.predictedTimeAveragedCoefficients = squeeze(mean(K.predictedCoefficients,2));
end

%% ------------------------------------------------------------------------
function fig = newFigure(name,opts)
fig = figure('Name',char(name),'NumberTitle','off', ...
  'Color','w','Visible',char(opts.Visible));
if opts.Docked, fig.WindowStyle = 'docked'; end
end

%% ------------------------------------------------------------------------
function plotFunctions(fig,G,allGroups,S,alloc,rawLimit,alignedLimit)
nAlloc=numel(alloc);
tl = tiledlayout(fig,3,nAlloc,'TileSpacing','compact','Padding','compact');
isCombined = G.name == "Combined";
members = unique(G.trials.animal,'stable');
for row = 1:2
  if row == 1
    field = 'totalStepCohPC';
    limit = rawLimit;
    fitField = 'rawFit';
    xlabelText = 'Total coherence step (%)';
  else
    field = 'alignedStep';
    limit = alignedLimit;
    fitField = 'alignedFit';
    xlabelText = 'Total step / session reference threshold';
  end
  for a = 1:nAlloc
    ax = nexttile(tl,(row-1)*nAlloc+a);
    hold(ax,'on');
    if isCombined
      for m = 1:numel(members)
        H = allGroups(find([allGroups.name] == members(m),1));
        drawBinned(ax,H.trials,field,alloc(a),limit,animalColor(members(m)),'o');
        drawFit(ax,H.(fitField),a,limit,animalColor(members(m)),'--',1.25);
      end
      drawFit(ax,G.(fitField),a,limit,[0 0 0],'-',2);
    else
      color = animalColor(G.name);
      drawBinned(ax,G.trials,field,alloc(a),limit,color,'o');
      drawFit(ax,G.(fitField),a,limit,color,'-',2);
    end
    yline(ax,.75,':','Color',[.25 .25 .25], ...
      'LineWidth',1.7,'HandleVisibility','off');
    if row == 2
      xline(ax,1,'-','Color',[.35 .35 .35], ...
        'LineWidth',1.5,'HandleVisibility','off');
    end
    xlim(ax,limit); ylim(ax,[.45 1.03]);
    grid(ax,'on');
    title(ax,sprintf('%g:%g',alloc(a),100-alloc(a)));
    if a == 1, ylabel(ax,'Fraction correct'); end
    xlabel(ax,xlabelText);
    if row == 1 && a == 1
      text(ax,.03,.05,'Raw steps','Units','normalized','FontWeight','bold');
      if isCombined
        text(ax,.03,.94,'Meetz blue; Neesha orange; pooled fit black', ...
          'Units','normalized','FontSize',8,'VerticalAlignment','top');
      end
    elseif row == 2 && a == 1
      text(ax,.03,.05,'Session aligned','Units','normalized','FontWeight','bold');
    end
  end
end
ax = nexttile(tl,2*nAlloc+1,[1 min(2,nAlloc)]);
hold(ax,'on');
R = [S.threshold100 S.threshold75 S.threshold62_5 S.threshold50];
[~,columns]=ismember(alloc,[100 75 62.5 50]); R=R(:,columns);
R = R ./ R(:,1);
yline(ax,1,':','Color',[.5 .5 .5],'HandleVisibility','off');
if isCombined
  for m = 1:numel(members)
    idx = S.animal == members(m);
    drawThresholdMean(ax,R(idx,:),animalColor(members(m)), ...
      char(members(m)),(m-1.5)*.1);
  end
  drawThresholdMean(ax,R(ismember(S.animal,members),:), ...
    [0 0 0],'Combined',0);
else
  idx = S.animal == G.name;
  drawThresholdMean(ax,R(idx,:),animalColor(G.name), ...
    char(G.name),0);
end
xticks(ax,1:nAlloc); xticklabels(ax,idlAllocationLabels(alloc));
xlim(ax,[.8 nAlloc+.2]); grid(ax,'on');
ylabel(ax,'75% threshold / same-session 100:0 threshold');
xlabel(ax,'Unsigned allocation');
legend(ax,'Location','northwest');
title(ax,'Session geometric mean and ±1 SE (log scale)');
sgtitle(tl,sprintf('%s | %d sessions, %d trials | %s', ...
  G.name+string(idlAllocationTag(G.twoAllocationMode)),numel(G.sessionIndex),height(G.trials),datestr(now,'yyyy-mm-dd')));
end

%% ------------------------------------------------------------------------
function drawThresholdMean(ax,ratios,color,label,offset)
logRatios = log(ratios);
meanLog = mean(logRatios,1);
meanRatio = exp(meanLog);
if size(ratios,1)>1
  seLog = std(logRatios,0,1)/sqrt(size(ratios,1));
else
  seLog = nan(1,size(ratios,2));
end
x = (1:size(ratios,2))+offset;
plot(ax,x,meanRatio,'o-','Color',color,'LineWidth',1.8, ...
  'MarkerFaceColor',color,'DisplayName',label);
lower = meanRatio(2:end)-exp(meanLog(2:end)-seLog(2:end));
upper = exp(meanLog(2:end)+seLog(2:end))-meanRatio(2:end);
errorbar(ax,x(2:end),meanRatio(2:end),lower,upper, ...
  'LineStyle','none','Color',color,'CapSize',5, ...
  'HandleVisibility','off');
end

%% ------------------------------------------------------------------------
function drawBinned(ax,T,field,alloc,limit,color,marker)
sel = abs(double(T.unsignedAllocPC)-alloc)<1e-6;
stepValues = double(T.(field));
x = stepValues(sel);
y = double(T.correct(sel));
edges = linspace(limit(1),limit(2),13);
bin = discretize(x,edges);
for b = 1:12
  index = bin == b;
  n = sum(index);
  if n < 8, continue; end
  success = sum(y(index));
  phat = success/n;
  z = 1.96;
  center = (phat+z^2/(2*n))/(1+z^2/n);
  half = z*sqrt(phat*(1-phat)/n+z^2/(4*n^2))/(1+z^2/n);
  errorbar(ax,mean(x(index)),phat,max(0,phat-center+half), ...
    max(0,center+half-phat),'Color',color,'LineStyle','none', ...
    'Marker',marker,'MarkerSize',4,'CapSize',3,'HandleVisibility','off');
end
end

%% ------------------------------------------------------------------------
function drawFit(ax,F,index,limit,color,style,width)
x = linspace(limit(1),limit(2),240);
prob = idqWeibullP(x,F.alpha(index),F.beta(index),F.lapse);
plot(ax,x,prob,'LineStyle',style,'Color',color, ...
  'LineWidth',width,'HandleVisibility','off');
end

%% ------------------------------------------------------------------------
function c = animalColor(name)
if name == "Meetz"
  c = [0.12 0.35 0.75];
elseif name == "Neesha"
  c = [0.85 0.35 0.14];
else
  c = [0.3 0.55 0.3];
end
end

%% ------------------------------------------------------------------------
function plotKernels(fig,G,signalFractionPC,scale,summaryScale)
nPanels=numel(signalFractionPC)+1; nColumns=min(4,nPanels);
tl = tiledlayout(fig,ceil(nPanels/nColumns),nColumns,'TileSpacing','compact','Padding','compact');
colors = [0.12 .35 .75;.85 .35 .14];
for a = 1:numel(signalFractionPC)
  ax = nexttile(tl,a);
  hold(ax,'on');
  preStep = G.kernel.preStepMS;
  xline(ax,0,'--','Color',[.35 .35 .35], ...
    'LineWidth',1.3,'HandleVisibility','off');
  yline(ax,0,'-','Color',[.72 .72 .72],'HandleVisibility','off');
  for side = 1:2
    kernel = G.kernel.coefficients(a,:,side);
    before = G.kernel.timeCentersMS<0;
    preMean = mean(kernel(before));
    stepMean = mean(kernel(~before));
    plot(ax,G.kernel.timeCentersMS,kernel,'-', ...
      'Color',colors(side,:),'LineWidth',1.8, ...
      'HandleVisibility','off');
    plot(ax,[-preStep 0],[preMean preMean],':', ...
      'Color',colors(side,:),'LineWidth',1.5,'HandleVisibility','off');
    plot(ax,[0 250],[stepMean stepMean],':', ...
      'Color',colors(side,:),'LineWidth',1.5,'HandleVisibility','off');
  end
  xlim(ax,[-preStep 250]);
  ylim(ax,[-1.08*scale 1.08*scale]);
  xticks(ax,[-preStep -500 -250 0 125 250]);
  grid(ax,'on');
  title(ax,sprintf('%g%% signal (n=%d)', ...
    signalFractionPC(a),G.kernel.nTrials(a)));
  if a>=5, xlabel(ax,'Time from step onset (ms)'); end
  if a==1 || a==5
    ylabel(ax,'Correct − error noise (% coherence)');
  end
end
ax = nexttile(tl,nPanels);
hold(ax,'on');
yline(ax,0,'-','Color',[.72 .72 .72],'HandleVisibility','off');
for side = 1:2
  plot(ax,signalFractionPC,G.kernel.relativeGain(:,side), ...
    'o-','Color',colors(side,:),'LineWidth',2, ...
    'MarkerSize',5, ...
    'DisplayName',sprintf('Observed %s',G.kernel.patchSource{side}));
end
modelColors = [.42 .42 .42;.48 .24 .62;.10 .54 .31;.75 .20 .55];
modelStyles = {':','-.','--','-'};
for model = 1:numel(G.kernel.predictedModelNames)
  plot(ax,signalFractionPC, ...
    G.kernel.predictedRelativeGain(:,1,model), ...
    'LineStyle',modelStyles{model},'Color',modelColors(model,:), ...
    'LineWidth',2, ...
    'DisplayName',sprintf('Predicted Change %s', ...
    G.kernel.predictedModelNames(model)));
end
xlim(ax,[0 100]);
ylim(ax,[-1.08*summaryScale 1.08*summaryScale]);
xticks(ax,[0 25 50 75 100]);
grid(ax,'on');
xlabel(ax,'Signal allocation (%)');
ylabel(ax,'Integrated gain / observed 100% Change gain');
if G.kernel.referenceGainAvailable
  title(ax,'Kernel-weighted noise regression');
else
  title(ax,'Relative gain unavailable: no 100% reference');
end
legend(ax,'Location','southeast','FontSize',8);
sgtitle(tl,sprintf('%s | choice-noise kernels | %d sessions, %d trials | %s', ...
  G.name+string(idlAllocationTag(G.twoAllocationMode)),numel(G.sessionIndex),height(G.trials),datestr(now,'yyyy-mm-dd')));
end
