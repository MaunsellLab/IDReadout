function sessionHeader = makeSessionHeader(header, trialMeta, trials)
% makeSessionHeader  Build stable analysis metadata for MT readout sessions.
%
% sessionHeader = makeSessionHeader(header, trialMeta, trials)
%
% INPUTS
%   header     Converted Lablib header struct.
%
%   trialMeta  Session-level metadata derived during conversion. Expected
%              fields include:
%                animal
%                nProbeDirections
%                probeDirectionsDeg
%                probeTags
%                nNoiseTrials
%
%              For header-only updates, the existing sessionHeader can be
%              passed as trialMeta because it contains these same fields.
%
%   trials     Parent-session trial cell array. When supplied, a compact
%              behaviorTrialTable is created with one row per parent trial.
%
% OUTPUT
%   sessionHeader
%       Session-level metadata plus:
%
%       behaviorTrialTable
%           One row per parent-session trial. This table contains only the
%           compact behavioral quantities needed for psychometric and gain
%           analyses; it does not contain stimulus time series or raw data.
%
%           Variables:
%             trialIdx       Parent-session trial index
%             trialEnd       Raw trialEnd code
%             trialCertify   Raw trialCertify code
%             isValid        trialCertify==0 and trialEnd is 0 or 1
%             correct        trialEnd==0 for valid trials; NaN otherwise
%             changeIndex    1=DEC, 2=INC
%             hasStepNoise   Logical cohNoise flag
%             baseCohPC      Baseline coherence (%)
%             stepCohPC      Step coherence (%)
%             stepDeltaCohPC stepCohPC - baseCohPC; 
%             prefCohNoisePC Preferred-stream noise amplitude (%)
%             probeCohNoisePC Probe-stream noise amplitude (%)
%             probeDirDeg    Probe offset (deg); NaN for no-noise trials
%             changeSide     0=RF, 1=Opp
%             chosenSide     0=RF, 1=Opp; NaN when unavailable
%
% The table preserves every parent trial so trialIdx remains an exact join
% key to derived probeSession files.

if nargin < 2 || isempty(trialMeta)
  trialMeta = struct();
end
if nargin < 3
  trials = [];
end

sessionHeader = struct();

% Parent-session probe-direction provenance. These describe the recording
% session as acquired, before trial correction, exclusion, or probe splitting.
sessionHeader.animal = getMetaField(trialMeta, 'animal', []);
sessionHeader.nProbeDirections = getMetaField(trialMeta, 'nProbeDirections', NaN);
sessionHeader.probeDirectionsDeg = getMetaField(trialMeta, 'probeDirectionsDeg', []);
sessionHeader.probeTags = getMetaField(trialMeta, 'probeTags', {});
sessionHeader.nNoiseTrials = getMetaField(trialMeta, 'nNoiseTrials', NaN);

% Stable fields used by downstream analyses or provenance displays.
copyFields = { ...
  'date', ...
  'subject', ...
  'taskName', ...
  'numberOfTrials', ...
  'frameRateHz', ...
  'preStepMS', ...
  'stepMS', ...
  'cohNoiseFrameMS', ...
  'prefDirDeg', ...
  'prefCohNoisePC' ...
  };

for k = 1:numel(copyFields)
  f = copyFields{k};
  if isfield(header, f)
    sessionHeader.(f) = localDataValue(header.(f));
  end
end

[~, sessionHeader.fileName] = fileparts(header.fileName);

% Modern files often store coherence-noise quantities in blockStatus.
if ~isfield(sessionHeader, 'prefCohNoisePC')
  [v, ok] = localGetBlockStatusValue(header, 'prefCohNoisePC');
  if ok
    sessionHeader.prefCohNoisePC = v;
  end
end
if ~isfield(sessionHeader, 'cohNoiseFrameMS')
  [v, ok] = localGetBlockStatusValue(header, 'cohNoiseFrameMS');
  if ok
    sessionHeader.cohNoiseFrameMS = v;
  end
end
if ~isfield(sessionHeader, 'prefCohNoisePC') && isfield(header, 'prefNoiseCohPC')
  sessionHeader.prefCohNoisePC = localDataValue(header.prefNoiseCohPC);
end

if ~isempty(trials)
  behaviorTrialTable = makeBehaviorTrialTable(trials);
  behaviorTrialTable.stepDeltaCohPC = behaviorTrialTable.stepCohPC - behaviorTrialTable.baseCohPC;
  sessionHeader.behaviorTrialTable = behaviorTrialTable;

  % These counts are useful consistency checks and make inventory creation
  % possible without loading the full trials cell array.
  T = sessionHeader.behaviorTrialTable;
  sessionHeader.nValidBehaviorTrials = sum(T.isValid);
  sessionHeader.nValidNoiseTrials = sum(T.isValid & T.hasStepNoise);
  sessionHeader.nValidNoNoiseTrials = sum(T.isValid & ~T.hasStepNoise);
end

sessionHeader.createdBy = mfilename;
sessionHeader.createdDate = datetime('now');
sessionHeader.metadataContract = [ ...
  'session-level analysis metadata plus compact behaviorTrialTable; ' ...
  'no raw trial data or stimulus time series'];
end

%% ------------------------------------------------------------------------
function T = makeBehaviorTrialTable(trials)
% makeBehaviorTrialTable  Extract compact, explicitly defined trial fields.

nTrials = numel(trials);

trialIdx       = (1:nTrials)';
trialEnd       = nan(nTrials, 1);
trialCertify   = nan(nTrials, 1);
isValid        = false(nTrials, 1);
correct        = nan(nTrials, 1);
changeIndex    = nan(nTrials, 1);
hasStepNoise   = false(nTrials, 1);
baseCohPC      = nan(nTrials, 1);
stepCohPC      = nan(nTrials, 1);
prefCohNoisePC = nan(nTrials, 1);
probeCohNoisePC = nan(nTrials, 1);
probeDirDeg    = nan(nTrials, 1);
changeSide     = nan(nTrials, 1);
chosenSide     = nan(nTrials, 1);

for iTrial = 1:nTrials
  tr = trials{iTrial};

  % Keep every parent trial. Missing or incomplete trial records remain
  % present with NaNs and isValid=false.
  if ~isstruct(tr)
    continue;
  end

  if isfield(tr, 'trialEnd') && isfield(tr.trialEnd, 'data')
    trialEnd(iTrial) = double(tr.trialEnd.data);
  end
  if isfield(tr, 'trialCertify') && isfield(tr.trialCertify, 'data')
    trialCertify(iTrial) = double(tr.trialCertify.data);
  end

  if isfield(tr, 'trial') && isfield(tr.trial, 'data')
    D = tr.trial.data;

    if isfield(D, 'changeIndex')
      changeIndex(iTrial) = double(D.changeIndex) + 1;  % 1=DEC, 2=INC
    end
    if isfield(D, 'cohNoise')
      hasStepNoise(iTrial) = logical(D.cohNoise);
    end
    if isfield(D, 'baseCohPC')
      baseCohPC(iTrial) = double(D.baseCohPC);
    end
    if isfield(D, 'stepCohPC')
      stepCohPC(iTrial) = double(D.stepCohPC);
    end
    if isfield(D, 'prefCohNoisePC')
      prefCohNoisePC(iTrial) = double(D.prefCohNoisePC);
    end
    if isfield(D, 'probeCohNoisePC')
      probeCohNoisePC(iTrial) = double(D.probeCohNoisePC);
    end
    if isfield(D, 'probeDirDeg')
      probeDirDeg(iTrial) = double(D.probeDirDeg);
    end
    if isfield(D, 'changeSide')
      changeSide(iTrial) = double(D.changeSide);
    end
  end


  if isfield(tr, 'targetChosen') && isfield(tr.targetChosen, 'data')
    chosenSide(iTrial) = double(tr.targetChosen.data);
  end

  isValid(iTrial) = ...
    trialCertify(iTrial) == 0 && ismember(trialEnd(iTrial), [0 1]);

  if isValid(iTrial)
    correct(iTrial) = double(trialEnd(iTrial) == 0);
  end
end

T = table( ...
  trialIdx, trialEnd, trialCertify, isValid, correct, changeIndex, hasStepNoise, baseCohPC, stepCohPC, prefCohNoisePC, ...
  probeCohNoisePC, probeDirDeg, changeSide, chosenSide);
end

%% ------------------------------------------------------------------------
function [v, ok] = localGetBlockStatusValue(header, fieldName)
ok = false;
v = [];
if isfield(header, 'blockStatus') && isfield(header.blockStatus, 'data') && ...
    isfield(header.blockStatus.data, fieldName)
  x = header.blockStatus.data.(fieldName);
  v = localExtractValue(x);
  ok = true;
end
end

%% ------------------------------------------------------------------------
function v = localDataValue(x)
v = x;
while isstruct(v) && isfield(v, 'data')
  v = v.data;
end
if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

%% ------------------------------------------------------------------------
function v = localExtractValue(x)
v = x;
while isstruct(v) && isfield(v, 'data')
  v = v.data;
end
if isnumeric(v) && ~isscalar(v)
  v = v(1);
end
end

%% ------------------------------------------------------------------------
function v = getMetaField(S, fieldName, defaultValue)
if isstruct(S) && isfield(S, fieldName)
  v = S.(fieldName);
else
  v = defaultValue;
end

if isnumeric(v)
  v = v(:)';
end
end
