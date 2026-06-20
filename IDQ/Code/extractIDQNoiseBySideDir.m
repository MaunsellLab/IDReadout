function noiseBySideDir = extractIDQNoiseBySideDir(trials, trialData, sessionHeader)
% extractIDQNoiseBySideDir
%
% Extract full six-stream IDQ noise for valid completed trials.
%
% Output:
%   noiseBySideDir: 2 sides x 3 directions x nFrames x nTrials
%
% Side with coherence increment step:
%   1 = RF, 2 = Opp
% Direction dimension:
%   1, 2, 3 corresponding to sessionHeader.dirDeg
%
% Noise is kept in raw task sign convention. DEC trials are not flipped.


msPerVFrame = 1000.0 / sessionHeader.frameRateHz;
m = round((sessionHeader.preStepMS + sessionHeader.stepMS) / msPerVFrame);
nTrials = numel(trialData.validIdx);
noiseBySideDir = nan(2, 3, m, nTrials);
t = trials(trialData.validIdx);

for kk = 1:nTrials
  tr = t{kk};

  changeNoise = fillFromTimes('change', tr, m, msPerVFrame);
  noChangeNoise = fillFromTimes('noChange', tr, m, msPerVFrame);

  if tr.trial.data.changeSide == 0
    % RF changed, Opp noChange
    noiseBySideDir(1, :, :, kk) = changeNoise;
    noiseBySideDir(2, :, :, kk) = noChangeNoise;
  else
    % Opp changed, RF noChange
    noiseBySideDir(1, :, :, kk) = noChangeNoise;
    noiseBySideDir(2, :, :, kk) = changeNoise;
  end
end
end

%%=========================================================================
function noise = fillFromTimes(fieldName, tr, m, msPerVFrame)

noise0 = tr.([fieldName, 'CohsPC0']).data;
noise1 = tr.([fieldName, 'CohsPC1']).data;
noise2 = tr.([fieldName, 'CohsPC2']).data;
timesMS = tr.([fieldName, 'TimesMS']).data;

nTimes = min(numel(timesMS), numel(noise0));
timesMS = timesMS(1:nTimes);
noise0  = noise0(1:nTimes);
noise1  = noise1(1:nTimes);
noise2  = noise2(1:nTimes);

noise = nan(1, 3, m);

for tIndex = 1:nTimes
  t0 = timesMS(tIndex);
  theVFrame = max(1, floor(t0 / msPerVFrame) + 1);
  if theVFrame > m
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
  noise(:, 1, theVFrame:nextVFrame - 1) = noise0(tIndex);
  noise(:, 2, theVFrame:nextVFrame - 1) = noise1(tIndex);
  noise(:, 3, theVFrame:nextVFrame - 1) = noise2(tIndex);
end
end
