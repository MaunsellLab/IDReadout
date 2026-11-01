function validateIDJTopUpExclusion(sourcePath)
% validateIDJTopUpExclusion  Inject flagged trials and verify early exclusion.
% The source file is read only. The altered FullSession lives in a temp folder.
base = makeIDJSessionAnalysis(sourcePath);
assert(idjHasTopUpExclusion(base));
assert(all(~base.trialTable.topUpTrial));
assert(height(base.trialTable)>=3);
raw = load(sourcePath,'trials','sessionHeader','header');
selected = unique([1 2 height(base.trialTable)]);
sourceIndices = base.trialTable.sourceTrialIndex(selected);
for k = sourceIndices(:)'
  raw.trials{k}.trial.data.topUpTrial = 1;
  % A flagged trial must be removed BEFORE step/configuration validation.
  raw.trials{k}.trial.data.totalStepCohPC = -9999;
end
folder = tempname;
mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
path = fullfile(folder,'flagged_FullSession.mat');
save(path,'-struct','raw','-v7.3');
filtered = makeIDJSessionAnalysis(path);
keep = ~ismember(base.trialTable.sourceTrialIndex,sourceIndices);
assert(idjHasTopUpExclusion(filtered));
assert(isequal(filtered.trialTable.sourceTrialIndex,base.trialTable.sourceTrialIndex(keep)));
assert(~any(ismember(filtered.streamTable.sourceTrialIndex,sourceIndices)));
assert(height(filtered.streamTable)==3*sum(keep));
assert(filtered.nTrialsTotal==base.nTrialsTotal);
assert(filtered.nAnalysisTrials==base.nAnalysisTrials-numel(selected));
assert(filtered.nExcludedTrials==base.nExcludedTrials+numel(selected));
assert(filtered.nExcludedTopUpTrials==base.nExcludedTopUpTrials+numel(selected));
assert(filtered.nCorrect==sum(base.trialTable.correct(keep)));
assert(sum(filtered.configurationTable.nAnalysisTrials)==sum(keep));
assert(sum(filtered.stepTable.nAnalysisTrials)==sum(keep));
fields = {'changeNoiseByFrameTrial','noChangeNoiseByFrameTrial', ...
  'fullChangeNoiseByFrameTrial','fullNoChangeNoiseByFrameTrial'};
for k = 1:numel(fields)
  expected = base.(fields{k});
  assert(isequal(filtered.(fields{k}),expected(:,keep,:)));
end
old = base;
old.version = 1;
assert(~idjHasTopUpExclusion(old));
old = base;
old.trialTable.topUpTrial(1) = true;
assert(~idjHasTopUpExclusion(old));
rejected = false;
try
  fitIDJSessionPsychometrics(old.trialTable);
catch err
  rejected = strcmp(err.identifier,'fitIDJSessionPsychometrics:TopUpTrials');
end
assert(rejected,'The standalone fitter must reject flagged input.');
fprintf('IDJ top-up exclusion passed: %d injected trials removed with aligned noise.\n', ...
  numel(selected));
end
