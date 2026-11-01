function staleProbeOffsets = makeIDSProbeSessions(varargin)
% makeIDSProbeSessions  Split IDS parent sessions by unsigned probe offset.
%
% Unlike the IDR splitter, this function retains noisy and no-noise trials.
% Instruction and invalid trials remain in probeTrials for provenance but are
% excluded later when analysis rows are constructed.

p = inputParser;
p.addParameter('Replace', false, @islogical);
p.addParameter('Verbose', true, @islogical);
p.addParameter('IncludeIDQ', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

idsRoot = domainFolder(mfilename('fullpath'));
sources = makeSources(idsRoot, opts.IncludeIDQ);
staleProbeOffsets = [];

for iSource = 1:numel(sources)
  source = sources(iSource);
  if ~isfolder(source.fullFolder)
    if opts.Verbose
      fprintf('makeIDSProbeSessions: skipping missing %s FullSessions folder %s\n', ...
        source.domain, source.fullFolder);
    end
    continue
  end
  dataFilePaths = selectAnalysisFiles(source.fullFolder, 'IncludeIDQ', false);

  for iFile = 1:numel(dataFilePaths)
    sourcePath = dataFilePaths{iFile};
    [~, sourceBase] = fileparts(sourcePath);
    S = load(sourcePath, 'header', 'sessionHeader', 'trials');
    if ~isfield(S, 'trials') || isempty(S.trials)
      continue
    end
    if ~isfield(S, 'sessionHeader')
      S.sessionHeader = makeParentHeader(S.header, sourceBase, numel(S.trials));
    end
    if strcmp(source.domain, 'IDQ')
      [S.header, S.sessionHeader, S.trials] = ...
        normalizeIDQFullSessionForIDS(S.header, S.sessionHeader, S.trials);
    end

    [trialOffsets, hasTrialDesc] = extractTrialOffsets(S.trials);
    offsets = unique(trialOffsets(hasTrialDesc & isfinite(trialOffsets)));
    if isempty(offsets)
      error('makeIDSProbeSessions:NoProbeOffsets', ...
        'No trial.data.probeOffsetDeg values found in %s.', sourceBase);
    end
    if strcmp(source.domain, 'IDQ')
      assert(isequal(offsets(:)', 120), ...
        'makeIDSProbeSessions:UnexpectedIDQProbeOffset', ...
        'Expected only the 120-degree IDQ probe offset in %s; found %s.', ...
        sourceBase, mat2str(offsets));
    end

    for iOffset = 1:numel(offsets)
      probeOffsetDeg = double(offsets(iOffset));
      validateattributes(probeOffsetDeg, {'numeric'}, {'scalar','integer','positive','<=',180});
      probeTag = sprintf('Probe%d', round(probeOffsetDeg));
      idx = find(hasTrialDesc & trialOffsets == probeOffsetDeg);
      probeTrials = S.trials(idx);
      trialIdx = idx(:)';

      [probeFolder, outputName] = probeSessionDestination(source, sourceBase, probeTag);
      probeFolder = validFolder(probeFolder);
      outPath = fullfile(probeFolder, outputName);
      if isfile(outPath) && ~opts.Replace
        continue
      end

      sessionHeader = S.sessionHeader;
      sessionHeader.sourceDomain = source.domain;
      sessionProbeHeader = makeProbeHeader(sourceBase, sourcePath, source.domain, ...
        probeTag, probeOffsetDeg, offsets, probeTrials, trialIdx);
      sessionProbeHeader.probeSessionPath = outPath;
      sessionProbeHeader.trialIdxContract = ...
        'trialIdx indexes original parent-session trials represented by probeTrials';
      header = S.header; %#ok<NASGU>
      save(outPath, 'header', 'sessionHeader', 'sessionProbeHeader', ...
        'probeTrials', 'trialIdx', '-v7.3');
      staleProbeOffsets(end + 1) = probeOffsetDeg; %#ok<AGROW>
      if opts.Verbose
        fprintf('makeIDSProbeSessions: saved %s (%d candidate trials)\n', ...
          outPath, numel(probeTrials));
      end
    end
  end
end

staleProbeOffsets = unique(staleProbeOffsets);
end

function sources = makeSources(idsRoot, includeIDQ)
sources = struct('domain', 'IDS', ...
  'root', idsRoot, ...
  'fullFolder', fullfile(idsRoot, 'Data', 'FullSessions'));
if includeIDQ
  domainsFolder = fileparts(idsRoot);
  idqRoot = fullfile(domainsFolder, 'IDQ');
  sources(end + 1) = struct('domain', 'IDQ', ...
    'root', idqRoot, ...
    'fullFolder', fullfile(idqRoot, 'Data', 'FullSessions'));
end
end

function [folder, fileName] = probeSessionDestination(source, sourceBase, probeTag)
if strcmp(source.domain, 'IDQ')
  folder = fullfile(source.root, 'Data', probeTag);
  fileName = sprintf('ProbeSession%s_%s.mat', sourceBase, probeTag);
else
  folder = fullfile(source.root, 'Data', probeTag, 'ProbeSessions');
  fileName = sprintf('%s_%s.mat', sourceBase, probeTag);
end
end

function [offsets, hasTrialDesc] = extractTrialOffsets(trials)
n = numel(trials);
offsets = nan(1, n);
hasTrialDesc = false(1, n);
for k = 1:n
  tr = trials{k};
  if isfield(tr, 'trial') && isfield(tr.trial, 'data') && ...
      isfield(tr.trial.data, 'probeOffsetDeg')
    hasTrialDesc(k) = true;
    offsets(k) = double(tr.trial.data.probeOffsetDeg);
  end
end
end

function H = makeProbeHeader(sourceBase, sourcePath, sourceDomain, probeTag, probeOffsetDeg, ...
    parentOffsets, probeTrials, trialIdx)
H = struct();
H.version = 1;
H.sessionID = sourceBase;
H.sourceSessionID = sourceBase;
H.sourceDomain = sourceDomain;
H.parentFilePath = sourcePath;
H.probeTag = probeTag;
H.probeOffsetDeg = probeOffsetDeg;
H.parentProbeOffsetsDeg = parentOffsets(:)';
H.parentNProbeOffsets = numel(parentOffsets);
H.parentTrialIdx = trialIdx(:)';
H.nCandidateTrials = numel(probeTrials);

n = numel(probeTrials);
isInstruction = false(1, n);
isValidEOT = false(1, n);
isCertified = false(1, n);
amplitude = nan(1, n);
nStreams = nan(1, n);
randomized = false(1, n);
randomBaseDir = nan(1, n);
for k = 1:n
  tr = probeTrials{k};
  D = tr.trial.data;
  isInstruction(k) = logical(D.instructTrial);
  isValidEOT(k) = isfield(tr, 'extendedEOT') && ismember(double(tr.extendedEOT.data), [0 1]);
  isCertified(k) = isfield(tr, 'trialCertify') && double(tr.trialCertify.data) == 0;
  amplitude(k) = double(D.cohNoisePC);
  randomized(k) = logical(D.randomizeDriftDir);
  if isfield(D, 'randomBaseDir')
    randomBaseDir(k) = logical(D.randomBaseDir);
  end
  if isfield(tr, 'blockStatus') && isfield(tr.blockStatus, 'data')
    nStreams(k) = double(tr.blockStatus.data(1).numNoiseStreams);
  end
end
analysisMask = ~isInstruction & isValidEOT & isCertified;
H.nInstructionTrials = sum(isInstruction);
H.nAnalysisTrials = sum(analysisMask);
H.nInvalidTrials = n - H.nAnalysisTrials - H.nInstructionTrials;
H.nNoiseTrials = sum(analysisMask & arrayfun(@(q) logical(probeTrials{q}.trial.data.stepHasNoise), 1:n));
H.nNoNoiseTrials = H.nAnalysisTrials - H.nNoiseTrials;
H.noiseAmplitudePC = uniqueScalar(amplitude(analysisMask), 'noise amplitude');
H.numNoiseStreams = uniqueScalar(nStreams(analysisMask & isfinite(nStreams)), 'number of noise streams');
H.randomizeDriftDir = all(randomized(analysisMask));
if strcmp(sourceDomain, 'IDQ')
  idqBaseFlags = randomBaseDir(analysisMask & isfinite(randomBaseDir));
  assert(~isempty(idqBaseFlags) && all(idqBaseFlags == idqBaseFlags(1)), ...
    'makeIDSProbeSessions:MixedIDQBaseDirectionMode', ...
    'IDQ session %s mixes fixed and randomized base-direction trials.', sourceBase);
  H.randomBaseDir = logical(idqBaseFlags(1));
  if H.randomBaseDir
    H.baseDirectionMode = 'randomized orientation; directions spaced by 120 degrees';
  else
    H.baseDirectionMode = 'fixed orientation; drift selected by dirIndex';
  end
end

if H.numNoiseStreams == 0
  H.numNoiseStreams = 3;
end
assert(H.numNoiseStreams == 3, 'makeIDSProbeSessions:ExpectedThreeStreams', ...
  'Expected exactly three noise streams in %s.', sourceBase);

assert(H.randomizeDriftDir, 'makeIDSProbeSessions:DriftNotRandomized', ...
  'Not every analyzable trial randomized drift direction in %s.', sourceBase);
H.streamNames = {'Drift','Drift+Probe','Drift-Probe'};
H.streamDirectionDefinition = ...
  '[driftDirDeg, mod(driftDirDeg+probeOffsetDeg,360), mod(driftDirDeg-probeOffsetDeg,360)]';
H.createdBy = mfilename;
H.createdAt = datetime('now');
end

function value = uniqueScalar(values, label)
values = unique(round(double(values(isfinite(values))), 6));
assert(isscalar(values), 'makeIDSProbeSessions:MixedMetadata', ...
  'Expected one %s value; found %s.', label, mat2str(values));
value = values;
end

function H = makeParentHeader(header, fileName, nTrials)
H = struct();
[~, H.fileName] = fileparts(fileName);
H.numberOfTrials = nTrials;
for name = {'date','frameRateHz','preStepMS','stepMS'}
  f = name{1};
  if isfield(header, f)
    H.(f) = localDataValue(header.(f));
  end
end
end

function value = localDataValue(value)
while isstruct(value) && isfield(value, 'data')
  value = value.data;
end
if isnumeric(value) && ~isscalar(value)
  value = value(1);
end
end
