function trialData = extractIDQTrialData(trials)
% extractIDQTrialData
%
% Extract analysis-facing trialData from raw converted IDQ trials.
%
% Raw trials may use task-native 0-based indexing.
% trialData uses MATLAB/IDR-facing 1-based indexing.
%
% Conventions:
%   trialData.correct         logical
%   trialData.stepSignIndex   1 = DEC, 2 = INC
%   trialData.sideIndex       1 = RF, 2 = Opp
%   trialData.chosenSideIndex 1 = RF, 2 = Opp
%   trialData.dirIndex        1, 2, 3

eotCode = getTrialField(trials, 'extendedEOT', 'data');
certify = getTrialField(trials, 'trialCertify', 'data');
validMask = (eotCode == 0 | eotCode == 1) & (certify == 0);
trialData.validIdx = find(validMask);
t = trials(trialData.validIdx);
trialData.correct = eotCode(trialData.validIdx) == 0;
trialData.trialCertify = getTrialField(t, 'trialCertify', 'data');
trialData.chosenSideIndex = getTrialField(t, 'targetChosen', 'data') + 1;
trialData.stepSignIndex = getTrialField(t, 'trial', 'data', 'stepIndex') + 1;
trialData.sideIndex = getTrialField(t, 'trial', 'data', 'changeSide') + 1;
trialData.dirIndex = getTrialField(t, 'trial', 'data', 'dirIndex') + 1;
trialData.stepCoh = getTrialField(t, 'trial', 'data', 'stepCohPC');
trialData.hasStepNoise = logical(getTrialField(t, 'trial', 'data', 'stepHasNoise'));
trialData.sideBias = getTrialField(t, 'RFBias', 'data');
trialData.dirBias = getTrialField(t, 'dirBias', 'data');

end