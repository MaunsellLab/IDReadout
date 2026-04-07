function [prefNoise, probeNoise, trialOutcomes, changeSides, changeIndices] = ...
  extractPatchNoiseMatrices(header, trials, stepTypes)
% extractPatchNoiseMatrices
% Return patchwise noise matrices for all valid trials.
%
% OUTPUTS
%   prefNoise      : 2 x m x nTrials   (1=RF, 2=Opp)
%   probeNoise     : 2 x m x nTrials   (1=RF, 2=Opp)
%   trialOutcomes  : 1 x nTrials       (0=correct, 1=wrong)
%   changeSides    : 1 x nTrials       (0=RF changed, 1=Opp changed)
%   changeIndices  : 1 x nTrials       (1=DEC, 2=INC)

nTrials = numel(trials);
if nTrials == 0
  error('extractPatchNoiseMatrices:EmptyInput', 'Input "trials" is empty.');
end

validIdx = [];
trialOutcomes = [];
changeSides = [];
changeIndices = [];

for k = 1:nTrials
  tr = trials{k};
  if ~isfield(tr, 'trialEnd') || ~isfield(tr, 'trialCertify') || ~isfield(tr, 'trial')
    continue;
  end

  tCert = tr.trialCertify.data;
  tEnd  = tr.trialEnd.data;
  tStep = tr.trial.data.changeIndex + 1;   % 1=DEC, 2=INC
  if ~(tCert == 0 && ismember(tEnd, [0 1]) && ismember(tStep, stepTypes))
    continue;
  end

  hasAnyNoise = ...
    ~(isNoNoise(tr.changePrefCohsPC.data(:))   && isNoNoise(tr.changeProbeCohsPC.data(:)) && ...
    isNoNoise(tr.noChangePrefCohsPC.data(:)) && isNoNoise(tr.noChangeProbeCohsPC.data(:)));
  if ~hasAnyNoise
    continue;
  end

  validIdx(end+1)       = k; %#ok<AGROW>
  trialOutcomes(end+1)  = tEnd; %#ok<AGROW>
  changeSides(end+1)    = tr.trial.data.changeSide; %#ok<AGROW>
  changeIndices(end+1)  = tStep; %#ok<AGROW>
end

nValid = numel(validIdx);
if nValid == 0
  error('extractPatchNoiseMatrices:NoValidTrials', 'No valid matching trials were found.');
end

msPerVFrame = 1000.0 / header.frameRateHz.data(1);
m = round((header.preStepMS.data(1) + header.stepMS.data(1)) / msPerVFrame);

prefNoise  = nan(2, m, nValid);
probeNoise = nan(2, m, nValid);

for kk = 1:nValid
  tr = trials{validIdx(kk)};

  prefChange    = fillFromTimes(tr.changePrefCohsPC.data(:),    tr.changeTimesMS.data(:),   m, msPerVFrame);
  probeChange   = fillFromTimes(tr.changeProbeCohsPC.data(:),   tr.changeTimesMS.data(:),   m, msPerVFrame);
  prefNoChange  = fillFromTimes(tr.noChangePrefCohsPC.data(:),  tr.noChangeTimesMS.data(:), m, msPerVFrame);
  probeNoChange = fillFromTimes(tr.noChangeProbeCohsPC.data(:), tr.noChangeTimesMS.data(:), m, msPerVFrame);

  if tr.trial.data.changeSide == 0
    % RF changed, Opp noChange
    prefNoise(1,:,kk)  = prefChange;
    probeNoise(1,:,kk) = probeChange;
    prefNoise(2,:,kk)  = prefNoChange;
    probeNoise(2,:,kk) = probeNoChange;
  else
    % Opp changed, RF noChange
    prefNoise(1,:,kk)  = prefNoChange;
    probeNoise(1,:,kk) = probeNoChange;
    prefNoise(2,:,kk)  = prefChange;
    probeNoise(2,:,kk) = probeChange;
  end
end
end

function v = fillFromTimes(cohsPC, timesMS, m, msPerVFrame)
v = nan(m,1);
if isempty(timesMS) || isempty(cohsPC)
  return;
end

nTimes = min(numel(timesMS), numel(cohsPC));
timesMS = timesMS(1:nTimes);
cohsPC  = cohsPC(1:nTimes);

for tIndex = 1:nTimes
  t0 = timesMS(tIndex);
  theVFrame = floor(t0 / msPerVFrame) + 1;
  if theVFrame < 1
    theVFrame = 1;
  elseif theVFrame > m
    continue;
  end

  if tIndex < nTimes
    t1 = timesMS(tIndex + 1);
    nextVFrame = floor(t1 / msPerVFrame) + 1;
  else
    nextVFrame = m + 1;
  end

  if nextVFrame <= theVFrame
    continue;
  end
  if nextVFrame > m + 1
    nextVFrame = m + 1;
  end

  v(theVFrame:nextVFrame-1) = cohsPC(tIndex);
end
end

function tf = isNoNoise(x)
tf = isscalar(x) && isequal(x, 0);
end