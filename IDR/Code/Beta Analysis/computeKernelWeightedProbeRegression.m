function reg = computeKernelWeightedProbeRegression(probeSessionPath, weightData)
% computeKernelWeightedProbeRegression  Fit preferred/probe readout betas.
%
% The same parent-session leave-one-out temporal weights are applied to the
% change-side preferred and effective yoked-probe noise streams. Positive
% predictors are signed to favor a correct response for both INC and DEC.

if nargin < 2 || isempty(weightData)
  baseFolder = domainFolder(mfilename('fullpath'));
  weightPath = fullfile(baseFolder, 'Data', 'FullSessions', ...
    'BetaAnalysis', 'AcrossSessions', 'BetaWeights.mat');
  S = load(weightPath, 'weightData');
  weightData = S.weightData;
else
  weightPath = '';
end

S = load(probeSessionPath, 'sessionHeader', 'sessionProbeHeader', ...
  'prefNoiseByPatch', 'probeNoiseByPatch', 'trialOutcomesAll', ...
  'changeSidesAll', 'changeIndicesAll');
required = {'sessionHeader','sessionProbeHeader','prefNoiseByPatch', ...
  'probeNoiseByPatch','trialOutcomesAll','changeSidesAll','changeIndicesAll'};
missing = required(~isfield(S, required));
if ~isempty(missing)
  error('computeKernelWeightedProbeRegression:MissingFields', ...
    '%s is missing: %s', probeSessionPath, strjoin(missing, ', '));
end

[parentWeightRow, parentFileName] = matchParentWeightRow( ...
  S.sessionHeader, S.sessionProbeHeader, weightData);
weights = double(weightData.leaveOneOutWeights(parentWeightRow, :));
stepTMS = double(weightData.stepTMS(:)');

[prefStep, stepFrameTMS] = extractStepFrames(S.prefNoiseByPatch, S.sessionHeader);
[probeStep, probeStepTMS] = extractStepFrames(S.probeNoiseByPatch, S.sessionHeader);
if numel(weights) ~= size(prefStep,2) || ...
    numel(weights) ~= size(probeStep,2)
  error('computeKernelWeightedProbeRegression:WeightLengthMismatch', ...
    'Weight count (%d) does not match step-frame count (%d).', ...
    numel(weights), size(prefStep,2));
end
if max(abs(stepFrameTMS - stepTMS)) > 1e-6 || ...
    max(abs(probeStepTMS - stepTMS)) > 1e-6
  error('computeKernelWeightedProbeRegression:WeightTimeMismatch', ...
    'Probe-session step-frame times do not match BetaWeights.stepTMS.');
end
if abs(sum(weights)-1) > 1e-12
  error('computeKernelWeightedProbeRegression:BadWeightNormalization', ...
    'Leave-one-session-out weights do not sum to one.');
end

prefByPatch = squeeze(sum(prefStep .* reshape(weights,1,[],1), 2));
probeByPatch = squeeze(sum(probeStep .* reshape(weights,1,[],1), 2));
if isvector(prefByPatch), prefByPatch = reshape(prefByPatch,2,[]); end
if isvector(probeByPatch), probeByPatch = reshape(probeByPatch,2,[]); end

changeSides = double(S.changeSidesAll(:));
changeIndex = double(S.changeIndicesAll(:));
correct = double(S.trialOutcomesAll(:) == 0);
changePatch = normalizePatchIndex(changeSides);
nTrials = numel(correct);

xPref = nan(nTrials,1);
xProbe = nan(nTrials,1);
for i = 1:nTrials
  xPref(i) = prefByPatch(changePatch(i), i);
  xProbe(i) = probeByPatch(changePatch(i), i);
end
stepSign = nan(nTrials,1);
stepSign(changeIndex == 1) = -1;
stepSign(changeIndex == 2) = 1;
if any(~isfinite(stepSign))
  error('computeKernelWeightedProbeRegression:BadChangeIndex', ...
    'changeIndicesAll must contain only 1=DEC or 2=INC.');
end
xPref = stepSign .* xPref;
xProbe = stepSign .* xProbe;

fitByStep = struct();
fitByStep.dec = fitSubset(changeIndex == 1, xPref, xProbe, correct, 'dec');
fitByStep.inc = fitSubset(changeIndex == 2, xPref, xProbe, correct, 'inc');
fitByStep.combined = fitSubset(true(nTrials,1), xPref, xProbe, correct, 'combined');

trialTable = table();
trialTable.correct = logical(correct);
trialTable.changeSide = changeSides;
trialTable.changeIndex = changeIndex;
trialTable.isInc = changeIndex == 2;
trialTable.stepSign = stepSign;
trialTable.effectivePrefNoisePC = xPref;
trialTable.effectiveProbeNoisePC = xProbe;

reg = struct();
reg.version = 3;
reg.analysisName = 'kernelWeightedProbeRegression';
reg.method = ['common parent-session leave-one-out preferred change-side ' ...
  'kernel weights applied to change-side preferred and effective yoked-probe noise'];
reg.probeSessionPath = probeSessionPath;
reg.sessionHeader = S.sessionHeader;
reg.sessionProbeHeader = S.sessionProbeHeader;
reg.parentSessionFileName = parentFileName;
reg.parentWeightRow = parentWeightRow;
reg.weightInfo = struct( ...
  'weightSourceFile', weightPath, ...
  'stepTMS', stepTMS, ...
  'weights', weights, ...
  'weightSum', sum(weights), ...
  'leaveOneOutKernelSum', weightData.leaveOneOutKernelSum(parentWeightRow));
reg.predictorInfo = struct( ...
  'units', 'percent coherence', ...
  'side', 'change', ...
  'signConvention', 'positive favors correct; INC=+1, DEC=-1', ...
  'probeConvention', ['probeNoiseByPatch is the effective yoked-pair ' ...
    'perturbation produced by makeProbeSessions']);
reg.trialTable = trialTable;
reg.fitByStep = fitByStep;
reg.nTrials = nTrials;
reg.nDec = sum(changeIndex == 1);
reg.nInc = sum(changeIndex == 2);
reg.createdBy = mfilename;
reg.createdDate = datetime('now');
end

function out = fitSubset(mask, xPref, xProbe, correct, label)
out = struct('label',label,'fitUsable',false,'message','');
try
  fit = fitFixedFloorNoiseRegression(xPref(mask), xProbe(mask), correct(mask));
  fit.label = label;
  out = fit;
catch ME
  out.label = label;
  out.fitUsable = false;
  out.message = ME.message;
  out.nTrials = sum(mask);
end
end

function [stepNoise, stepTMS] = extractStepFrames(noiseByPatch, sessionHeader)
if ndims(noiseByPatch) ~= 3 || size(noiseByPatch,1) ~= 2
  error('computeKernelWeightedProbeRegression:BadNoiseShape', ...
    'Noise matrix must be 2 x frames x trials.');
end
frameRateHz = headerScalar(sessionHeader, 'frameRateHz');
preStepMS = headerScalar(sessionHeader, 'preStepMS');
stepMS = headerScalar(sessionHeader, 'stepMS');
msPerFrame = 1000/frameRateHz;
nFrames = size(noiseByPatch,2);
tMS = (0:nFrames-1)*msPerFrame;
mask = tMS >= preStepMS & tMS < preStepMS + stepMS;
stepNoise = double(noiseByPatch(:,mask,:));
stepTMS = tMS(mask);
end

function idx = normalizePatchIndex(raw)
raw = raw(:);
vals = unique(raw(isfinite(raw)));
if all(ismember(vals,[0 1]))
  idx = raw + 1;
elseif all(ismember(vals,[1 2]))
  idx = raw;
else
  error('computeKernelWeightedProbeRegression:BadChangeSide', ...
    'changeSidesAll must contain 0/1 or 1/2 patch indices.');
end
end

function [row, parentName] = matchParentWeightRow(sessionHeader, probeHeader, W)
% candidates = strings(0,1);
% if isfield(probeHeader,'parentFileName')
%   candidates(end+1) = string(probeHeader.parentFileName); 
% end
% if isfield(sessionHeader,'fileName')
%   candidates(end+1) = string(sessionHeader.fileName); 
% end
% if isfield(probeHeader,'sessionID')
  % candidates(end+1) = string(probeHeader.sessionID); 
% end
weightNames = string(W.sessionFileNames(:));
% weightBases = strings(size(weightNames));
% for i = 1:numel(weightNames)
%   [~,b,e] = fileparts(weightNames(i));
%   weightBases(i) = b + e;
% end
matches = false(size(weightNames));
% for c = 1:numel(candidates)
%   % [~,b,e] = fileparts(candidates(c));
%   % fullCandidate = b + e;
%   % if strlength(e)==0, fullCandidate = b + ".mat"; end
%   % matches = matches | strcmpi(weightBases, fullCandidate);
matches = strcmpi(weightNames, probeHeader.sessionID);
% end
row = find(matches);
if numel(row) ~= 1
  error('computeKernelWeightedProbeRegression:ParentSessionMatch', ...
    'Could not uniquely match probe session to BetaWeights (%d matches).', numel(row));
end
% parentName = char(weightNames(row));

parentName = char(weightNames(row));
end

function v = headerScalar(H, fieldName)
if ~isfield(H,fieldName)
  error('computeKernelWeightedProbeRegression:MissingHeaderField', ...
    'sessionHeader.%s is required.', fieldName);
end
v = H.(fieldName);
if isstruct(v) && isfield(v,'data'), v = v.data; end
v = double(v(1));
end
