function M = makeIDRIDQMatchedSummary(varargin)
% makeIDRIDQMatchedSummary  Assemble compact, matched IDQ/IDR report data.
%
% This function loads already-fitted across-session summaries. It does not
% refit psychometrics, gains, or pooling models. The only new calculation is
% a descriptive IDQ physical-candidate overlap summary from the saved
% trialwise rectangular-step predictors.
%
% Name-value arguments:
%   IDQSummaryFolder     Folder containing the IDQ source MAT files
%   IDRSummaryFolder     Folder containing the IDR source MAT file
%   CommonSummaryFolder  Folder for the compact matched summary
%   SavePath       Compact output MAT path (default in CommonSummaryFolder)
%   SaveSummary    Save the compact structure (default true)

% Required source files:
%   IDQ_AcrossSideSummary.mat
%   IDQ_SignedNoiseSummary.mat
%   IDQ_NoChangeSignedNoiseSummary.mat
%   IDQ_OpponentPoolingSummary.mat
%   IDQ_OpponentPNormSummary.mat
%   IDQ_InteractionSummary.mat
%   IDR_SideGainAnalysis.mat


p = inputParser;
addParameter(p,'IDQSummaryFolder',defaultDomainSummaryFolder('IDQ'), ...
  @(x) ischar(x) || isstring(x));
addParameter(p,'IDRSummaryFolder',defaultDomainSummaryFolder('IDR'), ...
  @(x) ischar(x) || isstring(x));
addParameter(p,'CommonSummaryFolder',defaultDomainSummaryFolder('Common Code'), ...
  @(x) ischar(x) || isstring(x));
addParameter(p,'SavePath','',@(x) ischar(x) || isstring(x));
addParameter(p,'SaveSummary',true,@(x) islogical(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;

idqSummaryFolder = char(opts.IDQSummaryFolder);
idrSummaryFolder = char(opts.IDRSummaryFolder);
commonSummaryFolder = char(opts.CommonSummaryFolder);
if strlength(string(opts.SavePath)) == 0
  savePath = fullfile(commonSummaryFolder,'IDRIDQ_MatchedSummary.mat');
else
  savePath = char(opts.SavePath);
end

A = loadRequired(idqSummaryFolder,'IDQ_AcrossSideSummary.mat','acrossSummary');
SC = loadRequired(idqSummaryFolder,'IDQ_SignedNoiseSummary.mat','signedNoiseSummary');
SN = loadRequired(idqSummaryFolder,'IDQ_NoChangeSignedNoiseSummary.mat', ...
  'noChangeSignedNoiseSummary');
PO = loadRequired(idqSummaryFolder,'IDQ_OpponentPoolingSummary.mat', ...
  'opponentPoolingSummary');
PN = loadRequired(idqSummaryFolder,'IDQ_OpponentPNormSummary.mat', ...
  'opponentPNormSummary');
IX = loadRequired(idqSummaryFolder,'IDQ_InteractionSummary.mat', ...
  'interactionSummary');
R = loadRequired(idrSummaryFolder,'IDR_SideGainAnalysis.mat','inventory');

validateSources(A,SC,SN,PO,PN,IX,R);

M = struct();
M.createdAt = datetime('now');
M.createdBy = mfilename;
M.analysisType = 'matched descriptive IDQ-IDR report; no behavioral refitting';
M.idqSummaryFolder = string(idqSummaryFolder);
M.idrSummaryFolder = string(idrSummaryFolder);
M.commonSummaryFolder = string(commonSummaryFolder);
M.sourceFiles = string({ ...
  'IDQ_AcrossSideSummary.mat';'IDQ_SignedNoiseSummary.mat'; ...
  'IDQ_NoChangeSignedNoiseSummary.mat';'IDQ_OpponentPoolingSummary.mat'; ...
  'IDQ_OpponentPNormSummary.mat';'IDQ_InteractionSummary.mat'; ...
  'IDR_SideGainAnalysis.mat'});

% IDQ: retain fitted and binned quantities, but omit duplicated trial tables.
Q = struct();
Q.nSessions = A.nSessions;
Q.nTrials = A.nTrials;
Q.nStepNoiseTrials = A.nStepNoiseTrials;
Q.predictorType = string(A.noisePredictorType);
Q.psychometric = makeIDQPsychometricRecord(A);
Q.kernel.change = A.side.change.kernel;
Q.kernel.noChange = A.side.noChange.kernel;
Q.kernel.changeRelative = A.side.change.relKernels;
Q.kernel.noChangeRelative = A.side.noChange.relKernels;
Q.linear.change = A.side.change.noiseGain;
Q.linear.noChange = A.side.noChange.noiseGain;
Q.signed.change = removeFieldIfPresent(SC.signedFits,'trialData');
Q.signed.changeDriftCurve = SC.driftCurve;
Q.signed.changeNonDriftCurve = SC.nonDriftCurve;
Q.signed.noChange = removeFieldIfPresent(SN.signedFits,'trialData');
Q.signed.noChangeAllStreamCurve = SN.allStreamCurve;
Q.signed.noChangeRelativeCurves = SN.relativeCurves;
Q.pooling = removeFieldIfPresent(PO.poolingFits,'trialData');
Q.poolingDisagreementCurve = PO.disagreementCurve;
Q.pNorm = removeFieldIfPresent(PN.pNormFits,'trialData');
Q.interaction = removeFieldIfPresent(IX.interactionFits,'trialData');
Q.interactionSurface = IX.surfaceSummary;
Q.preferredNullRatio = PO.preferredToNullRatio;
Q.overlap = makeIDQPhysicalOverlap( ...
  SC.signedFits.trialData,SN.signedFits.trialData);
Q.stepNoiseSD = [std(SC.signedFits.trialData.nD,0,'omitnan'), ...
  std(SC.signedFits.trialData.nPlus,0,'omitnan'), ...
  std(SC.signedFits.trialData.nMinus,0,'omitnan')];
Q.meanStepCoh = mean(SC.signedFits.trialData.stepCoh,'omitnan');
M.idq = Q;

% IDR: tables are already compact and preserve offset-specific results.
D = struct();
D.nSessions = R.alignedPsychometric.primary.nSessions;
D.nTrials = R.alignedPsychometric.primary.nTrials;
D.nStepNoiseTrials = height(R.sideGainData.trialTable);
D.predictorType = "rectStep";
D.psychometric = R.alignedPsychometric.primary;
D.sessionPsychFits = R.sessionPsychFits;
D.commonKernel.tMS = R.sideGainPredictors.tMS;
D.commonKernel.kernel = R.sideGainPredictors.pooledCommonKernel;
D.commonKernel.rectReference = R.sideGainPredictors.rectStep.weights;
D.commonKernel.stepFrames = R.sideGainPredictors.stepFrames;
D.predictorDiagnostics = R.sideGainPredictors.predictorDiagnostics;
D.linearSide = R.sideGainFits.rectStep;
D.matchedGains = R.matchedGainFits.summaryTable;
D.signedGains = R.signedGainFits.summaryTable;
D.candidateOverlap = R.candidateOverlap.offsetSummary;
D.candidateOverlapBySession = R.candidateOverlap.sessionOffsetSummary;
D.pooling.offset = R.controlledPoolingFits.offsetSummary;
D.pooling.overall = R.controlledPoolingFits.overallSummary;
D.pooling.profile = R.controlledPoolingFits.sharedPNorm.profile.table;
D.pooling.sharedP = R.controlledPoolingFits.sharedPNorm.p;
D.pooling.profileLowerOneSided95 = ...
  R.controlledPoolingFits.sharedPNorm.profile.lowerOneSided95;
D.pooling.profileLowerTwoSided95 = ...
  R.controlledPoolingFits.sharedPNorm.profile.lowerTwoSided95;
D.preferredNullRatio = R.controlledPoolingFits.preferredNullRatio;
D.excludedProbeDirectionsDeg = R.sideGainData.excludedProbeDirectionsDeg;
D.probeScalingContract = string(R.sideGainData.probeScalingContract);
D.commonKernelDefinition = string(R.sideGainData.commonKernelDefinition);
M.idr = D;

M.analysisContract = makeAnalysisContract(M);

if opts.SaveSummary
  M.savePath = string(savePath);
  save(savePath,'M','-v7.3');
end
end

%% ------------------------------------------------------------------------
function folder = defaultDomainSummaryFolder(domainName)
try
  root = domainFolder(mfilename('fullpath'),domainName);
  folder = fullfile(root,'Data','AcrossSessionSummaries');
catch
  folder = pwd;
end
end

function value = loadRequired(folder,fileName,variableName)
path = fullfile(folder,fileName);
if ~isfile(path)
  error('makeIDRIDQMatchedSummary:MissingFile','Missing source file: %s',path);
end
S = load(path,variableName);
if ~isfield(S,variableName)
  error('makeIDRIDQMatchedSummary:MissingVariable', ...
    '%s does not contain %s.',path,variableName);
end
value = S.(variableName);
end

function validateSources(A,SC,SN,PO,PN,IX,R)
requireFields(A,{'psychFit','trialTable','side','noisePredictorType'},'IDQ across summary');
requireFields(SC,{'signedFits','driftCurve','nonDriftCurve'},'IDQ signed change');
requireFields(SN,{'signedFits','allStreamCurve'},'IDQ signed no-change');
requireFields(PO,{'poolingFits','preferredToNullRatio'},'IDQ pooling');
requireFields(PN,{'pNormFits'},'IDQ p-norm');
requireFields(IX,{'interactionFits','surfaceSummary'},'IDQ interactions');
requireFields(R,{'alignedPsychometric','sessionPsychFits','sideGainData', ...
  'sideGainPredictors','sideGainFits','matchedGainFits','signedGainFits', ...
  'candidateOverlap','controlledPoolingFits'},'IDR inventory');
if string(A.noisePredictorType) ~= "rectStep" || ...
    string(SC.noisePredictorType) ~= "rectStep" || ...
    string(SN.noisePredictorType) ~= "rectStep"
  error('makeIDRIDQMatchedSummary:PredictorMismatch', ...
    'Matched report requires rectangular-step IDQ predictors.');
end
end

function requireFields(S,names,label)
missing = names(~isfield(S,names));
if ~isempty(missing)
  error('makeIDRIDQMatchedSummary:MissingField', ...
    '%s lacks: %s.',label,strjoin(missing,', '));
end
end

function S = removeFieldIfPresent(S,name)
if isfield(S,name)
  S = rmfield(S,name);
end
end

%% ------------------------------------------------------------------------
function P = makeIDQPsychometricRecord(A)
T = A.trialTable;
x = double(T.alignedCoh(:));
y = double(T.correct(:));
valid = isfinite(x) & isfinite(y);
x = x(valid); y = y(valid);
edges = linspace(min(x),max(x),13);
edges(1) = -Inf; edges(end) = Inf;
bin = discretize(x,edges);
nBins = numel(edges)-1;
alignedCoh = nan(nBins,1); pCorrect = nan(nBins,1); nTrials = zeros(nBins,1);
for k = 1:nBins
  use = bin==k; nTrials(k)=sum(use);
  if any(use)
    alignedCoh(k)=mean(x(use),'omitnan');
    pCorrect(k)=mean(y(use),'omitnan');
  end
end
P = struct();
P.fit = A.psychFit.alignedWeibull;
P.binned = table(alignedCoh,pCorrect,nTrials);
P.sessionFits = A.psychFit.sessionFits;
P.nTrials = numel(y);
P.nSessions = A.nSessions;
ratio = double(T.stepCoh(:))./double(T.alignedCoh(:));
P.medianSessionThresholdCoh = median(ratio(isfinite(ratio)),'omitnan');
end

%% ------------------------------------------------------------------------
function O = makeIDQPhysicalOverlap(changeT,noChangeT)
required = {'nD','nPlus','nMinus'};
requireTableVariables(changeT,[required {'stepCoh'}],'IDQ change overlap');
requireTableVariables(noChangeT,required,'IDQ no-change overlap');

changePreferred = double(changeT.nD(:));
changeProbeMax = max(double(changeT.nPlus(:)),double(changeT.nMinus(:)));
step = double(changeT.stepCoh(:));
changeMargin = changeProbeMax-(step+changePreferred);

noChangePreferred = double(noChangeT.nD(:));
noChangeProbeMax = max(double(noChangeT.nPlus(:)),double(noChangeT.nMinus(:)));
noChangeMargin = noChangeProbeMax-noChangePreferred;

O = struct();
O.change = summarizeMargin(changeMargin,changePreferred,changeProbeMax,step);
O.noChange = summarizeMargin(noChangeMargin,noChangePreferred,noChangeProbeMax,[]);
O.change.histogram = makeHistogram(changeMargin,50);
O.noChange.histogram = makeHistogram(noChangeMargin,50);
O.changeDefinition = 'max(nPlus,nMinus) - (stepCoh + nD)';
O.noChangeDefinition = 'max(nPlus,nMinus) - nD';
O.note = ['Physical-coherence candidate overlap only; excludes MT tuning, ' ...
  'normalization, internal variability, and readout pooling.'];
end

function S = summarizeMargin(margin,preferred,probe,step)
valid = isfinite(margin) & isfinite(preferred) & isfinite(probe);
if ~isempty(step), valid = valid & isfinite(step); end
margin=margin(valid); preferred=preferred(valid); probe=probe(valid);
S = struct();
S.nTrials = numel(margin);
S.sdPreferred = std(preferred,0,'omitnan');
S.sdProbeMax = std(probe,0,'omitnan');
S.correlation = corr(preferred,probe,'Rows','complete');
S.meanMargin = mean(margin,'omitnan');
S.sdMargin = std(margin,0,'omitnan');
S.q50 = quantile(margin,0.50);
S.q90 = quantile(margin,0.90);
S.q95 = quantile(margin,0.95);
S.q99 = quantile(margin,0.99);
S.maximum = max(margin);
S.crossingProbability = mean(margin>0);
if isempty(step)
  S.meanStep = nan;
else
  S.meanStep = mean(step(valid),'omitnan');
end
end

function H = makeHistogram(x,nBins)
x=x(isfinite(x));
[count,edges]=histcounts(x,nBins,'Normalization','probability');
H.edges=edges(:); H.centers=0.5*(edges(1:end-1)+edges(2:end));
H.probability=count(:);
end

function requireTableVariables(T,names,label)
missing=setdiff(names,T.Properties.VariableNames);
if ~isempty(missing)
  error('makeIDRIDQMatchedSummary:MissingTableVariable', ...
    '%s lacks: %s.',label,strjoin(missing,', '));
end
end

%% ------------------------------------------------------------------------
function C = makeAnalysisContract(M)
C = struct();
C.common = { ...
  'INC trials'; ...
  'Physical coherence signs retained'; ...
  'Step-rectangular trialwise noise predictors'; ...
  'Session-threshold alignment'; ...
  'Pooled Weibull shape and lapse'; ...
  'Bernoulli trial likelihood'; ...
  'Winning mechanism remains behaviorally latent'};
C.idq = { ...
  sprintf('%d sessions; %d step-noise trials',M.idq.nSessions,M.idq.nStepNoiseTrials); ...
  'Three directions per patch: drift and +/-120 deg'; ...
  'Three physical noise streams per patch'; ...
  'Step magnitude is the physical coherence step'; ...
  sprintf('Step-integrated stream SDs: %.2f, %.2f, %.2f %% coherence',M.idq.stepNoiseSD); ...
  sprintf('Mean coherence step: %.2f %%',M.idq.meanStepCoh); ...
  sprintf('Preferred:null sensitivity ratio %.2f',M.idq.preferredNullRatio)};
C.idr = { ...
  sprintf('%d sessions; %d step-noise trials',M.idr.nSessions,M.idr.nStepNoiseTrials); ...
  'Preferred direction plus offset probes'; ...
  'Paired yoked probe candidates except 180 deg'; ...
  'Effective probe predictor sums yoked streams'; ...
  'Candidate predictor represents one probe mechanism'; ...
  'Step magnitude = stepCoh - prestepCoh'; ...
  sprintf('Excluded probe offset: %g deg',M.idr.excludedProbeDirectionsDeg); ...
  sprintf('Preferred:null sensitivity ratio %.2f',M.idr.preferredNullRatio)};
end
