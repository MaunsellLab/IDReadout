function [header, sessionHeader, trials] = normalizeIDQFullSessionForIDS(header, sessionHeader, trials)
% normalizeIDQFullSessionForIDS  Add the IDS FullSession contract to IDQ data.
%
% IDQ stores coherence traces by absolute stream index (0, 1, 2), whereas
% IDS stores them by their role on each trial (drift, drift+probe, and
% drift-probe). This function preserves the original IDQ fields and adds the
% IDS-compatible aliases and trial-descriptor fields needed by the IDS
% ProbeSession pipeline.

conversionVersion = 1;

if ~isfield(sessionHeader, 'animal') || isempty(sessionHeader.animal)
  animalNames = {'Meetz','Neesha'};
  matches = cellfun(@(s) contains(char(sessionHeader.fileName), s), animalNames);
  if sum(matches) == 1
    sessionHeader.animal = animalNames{matches};
  else
    sessionHeader.animal = 'unknown';
  end
end
sessionHeader.sourceDomain = 'IDQ';

if isfield(sessionHeader, 'idsFullSessionCompatibilityVersion')
  assert(sessionHeader.idsFullSessionCompatibilityVersion == conversionVersion, ...
    'normalizeIDQFullSessionForIDS:UnsupportedVersion', ...
    'Unsupported IDQ-to-IDS FullSession compatibility version %g.', ...
    double(sessionHeader.idsFullSessionCompatibilityVersion));
  validateNormalizedTrials(trials);
  return
end

randomBaseFlags = nan(1, numel(trials));
for t = 1:numel(trials)
  tr = trials{t};
  if ~isfield(tr, 'trial') || ~isfield(tr.trial, 'data')
    continue
  end
  D = tr.trial.data;
  required = {'changeSide','stepCohPC','dirIndex','dirDeg','cohNoisePC'};
  missing = required(~isfield(D, required));
  assert(isempty(missing), 'normalizeIDQFullSessionForIDS:MissingTrialField', ...
    'Trial %d lacks required IDQ field(s): %s.', t, strjoin(missing, ', '));

  dirIndex = double(D.dirIndex);
  assert(isscalar(dirIndex) && isfinite(dirIndex) && dirIndex == round(dirIndex) && ...
    ismember(dirIndex, 0:2), 'normalizeIDQFullSessionForIDS:InvalidDirIndex', ...
    'Trial %d has invalid dirIndex %s.', t, mat2str(dirIndex));

  if isfield(D, 'baseDirDeg')
    expectedDriftDir = mod(double(D.baseDirDeg) + 120 * dirIndex, 360);
    directionError = abs(mod(expectedDriftDir - double(D.dirDeg) + 180, 360) - 180);
    assert(directionError < 1e-6, ...
      'normalizeIDQFullSessionForIDS:DirectionMismatch', ...
      ['Trial %d has baseDirDeg %.6g and dirIndex %d, which predict drift ' ...
       'direction %.6g rather than dirDeg %.6g.'], ...
      t, double(D.baseDirDeg), dirIndex, expectedDriftDir, double(D.dirDeg));
  end

  % Preserve IDQ's increment explicitly before adopting IDS's endpoint
  % convention for stepCohPC.
  D.stepIncrementPC = D.stepCohPC;
  D.instructTrial = false;
  D.stepSign = 1;
  D.stepSide = D.changeSide;
  D.baseCohPC = D.cohNoisePC;
  D.stepCohPC = double(D.stepIncrementPC) + double(D.cohNoisePC);
  D.randomizeDriftDir = true;
  if ~isfield(D, 'randomBaseDir')
    D.randomBaseDir = false;
  end
  randomBaseFlags(t) = logical(D.randomBaseDir);
  D.driftDirDeg = D.dirDeg;
  D.probeOffsetDeg = 120;
  D.offsetIndex = 0;
  tr.trial.data = D;

  % IDQ stream indices increase in the positive direction around the three
  % directions. Thus dirIndex+1 is drift+120 and dirIndex+2 is drift-120.
  roleNames = {'Drift','Probe0','Probe1'};
  streamIndices = [dirIndex, mod(dirIndex + 1, 3), mod(dirIndex + 2, 3)];
  for side = {'change','noChange'}
    sideName = side{1};
    for r = 1:numel(roleNames)
      sourceField = sprintf('%sCohsPC%d', sideName, streamIndices(r));
      targetField = sprintf('%s%sCohsPC', sideName, roleNames{r});
      assert(isfield(tr, sourceField), ...
        'normalizeIDQFullSessionForIDS:MissingCoherenceTrace', ...
        'Trial %d lacks %s.', t, sourceField);
      tr.(targetField) = tr.(sourceField);
    end
  end

  if isfield(tr, 'blockStatus') && isfield(tr.blockStatus, 'data')
    B = tr.blockStatus.data;
    for b = 1:numel(B)
      B(b).randomizeDriftDir = true;
      B(b).numProbeOffsets = 1;
      B(b).probeDirDeg = 120;
    end
    tr.blockStatus.data = B;
  end
  trials{t} = tr;
end

sessionHeader.numProbeOffsets = 1;
sessionHeader.probeOffsetsDeg = 120;
sessionHeader.randomizeDriftDir = true;
baseModes = unique(randomBaseFlags(isfinite(randomBaseFlags)));
assert(isscalar(baseModes), ...
  'normalizeIDQFullSessionForIDS:MixedBaseDirectionMode', ...
  'An IDQ FullSession cannot mix fixed and randomized base-direction trials.');
sessionHeader.randomBaseDir = logical(baseModes);
if sessionHeader.randomBaseDir
  sessionHeader.baseDirectionMode = ...
    'randomized orientation; directions spaced by 120 degrees';
else
  sessionHeader.baseDirectionMode = ...
    'fixed orientation; drift selected by dirIndex';
end
sessionHeader.idsFullSessionCompatibilityVersion = conversionVersion;
sessionHeader.idsFullSessionCompatibility = ...
  ['Original IDQ fields retained; IDS role-based coherence aliases and ' ...
   'trial-descriptor fields added.'];

conversion = struct();
conversion.applied = true;
conversion.conversionVersion = conversionVersion;
conversion.appliedDate = datetime('now');
conversion.stepConvention = ...
  'stepIncrementPC retains IDQ increment; stepCohPC is pedestal plus increment';
conversion.streamConvention = ...
  'IDQ streams rotated by dirIndex into Drift, Probe0 (+120), Probe1 (-120)';
conversion.fixedOrientationConvention = ...
  'Absent randomBaseDir is interpreted as false';
if ~isfield(header, 'dataConversions') || isempty(header.dataConversions)
  header.dataConversions = struct();
end
header.dataConversions.IDQToIDSFullSession = conversion;

validateNormalizedTrials(trials);
end

function validateNormalizedTrials(trials)
requiredTrialFields = {'instructTrial','stepSign','stepSide','baseCohPC', ...
  'stepCohPC','stepIncrementPC','randomizeDriftDir','driftDirDeg', ...
  'probeOffsetDeg','offsetIndex'};
requiredTraceFields = {'changeDriftCohsPC','changeProbe0CohsPC', ...
  'changeProbe1CohsPC','noChangeDriftCohsPC','noChangeProbe0CohsPC', ...
  'noChangeProbe1CohsPC'};
for t = 1:numel(trials)
  tr = trials{t};
  if ~isfield(tr, 'trial') || ~isfield(tr.trial, 'data')
    continue
  end
  missingTrial = requiredTrialFields(~isfield(tr.trial.data, requiredTrialFields));
  missingTrace = requiredTraceFields(~isfield(tr, requiredTraceFields));
  assert(isempty(missingTrial) && isempty(missingTrace), ...
    'normalizeIDQFullSessionForIDS:IncompleteNormalizedTrial', ...
    'Normalized trial %d is missing required IDS-compatible fields.', t);
end
end
