function probeSessions = splitTrialsByProbeDirection(header, trials)
% splitTrialsByProbeDirection  Split one recording session into probe-specific analysis sessions.
%
% Each output element has:
%   .probeDirDeg
%   .probeTag
%   .header
%   .trials
%   .trialIdx
%
% Old single-probe files remain compatible if header.probeDirDeg.data exists
% and per-trial probeDirDeg is absent.

nTrials = numel(trials);
trialProbeDirs = nan(1, nTrials);

% ---- Preferred new format: selected probe direction stored per trial ----
for t = 1:nTrials
  if isfield(trials{t}, 'trial') && ...
      isfield(trials{t}.trial, 'data') && ...
      isfield(trials{t}.trial.data, 'probeDirDeg')

    trialProbeDirs(t) = double(trials{t}.trial.data.probeDirDeg);
  end
end

% ---- Backward compatibility: old single-probe sessions ----
if all(isnan(trialProbeDirs))
  if isfield(header, 'probeDirDeg') && isfield(header.probeDirDeg, 'data')
    probeDirDeg = double(header.probeDirDeg.data);
    trialProbeDirs(:) = probeDirDeg;
  else
    error('splitTrialsByProbeDirection:MissingProbeDir', ...
      ['No per-trial trials{t}.trial.data.probeDirDeg found, and no ' ...
       'header.probeDirDeg.data field is available for old-format compatibility.']);
  end
end

% ---- Assertions: fail early if partly encoded or malformed ----
if any(isnan(trialProbeDirs))
  bad = find(isnan(trialProbeDirs), 1, 'first');
  error('splitTrialsByProbeDirection:IncompleteProbeDir', ...
    'Missing trial.data.probeDirDeg for trial %d of %d.', bad, nTrials);
end
if any(~isfinite(trialProbeDirs))
  error('splitTrialsByProbeDirection:BadProbeDir', ...
    'Non-finite probe direction values found.');
end

% Exclude no-noise trials: probeDirDeg == -1 is the experimental-code
% sentinel for trials without coherence noise, not an analysis probe direction.
probeDirs = unique(trialProbeDirs);
probeDirs(probeDirs == -1) = [];
if isempty(probeDirs)
  error('splitTrialsByProbeDirection:NoProbeNoiseTrials', ...
    'No valid probe directions found after excluding probeDirDeg == -1 no-noise trials.');
end

probeSessions = repmat(struct( ...
  'probeDirDeg', [], ...
  'probeTag', '', ...
  'header', [], ...
  'trials', [], ...
  'trialIdx', []), 1, numel(probeDirs));

for p = 1:numel(probeDirs)
  probeDirDeg = probeDirs(p);
  idx = find(trialProbeDirs == probeDirDeg);

  probeHeader = header;
  probeHeader.probeDirDeg = struct('data', probeDirDeg);

  % Preserve the full randomized probe set if present, but make the selected
  % analysis-session probe direction authoritative via header.probeDirDeg.data.
  probeHeader.parentNumberOfTrials = nTrials;
  probeHeader = header;
  probeHeader.probeDirDeg = struct('data', probeDirDeg);

  probeSessions(p).probeDirDeg = probeDirDeg;
  probeSessions(p).probeTag = sprintf('probe%d', round(probeDirDeg));
  probeSessions(p).header = probeHeader;
  probeSessions(p).trials = trials(idx);
  probeSessions(p).trialIdx = idx;
end
end