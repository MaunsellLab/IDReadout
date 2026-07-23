function expected=fitIDRExpectedChoiceDOGGains(context,varargin)
% fitIDRExpectedChoiceDOGGains  Fit IDR gains to mean synthetic choices.
%
% Trial correctness is the fraction correct across internal-noise draws.
% Since the Bernoulli log likelihood is affine in the outcome, fitting these
% fractional outcomes equals fitting the mean log likelihood across draws.

context=validateAnalysisContext(context);
if context.Mode~="synthetic"
  error('IDReadout:SyntheticContextRequired', ...
    'Expected-choice fitting requires a synthetic context.');
end
validateSyntheticManifest(context);

p=inputParser;
addParameter(p,'SimulationPath',defaultSimulationPath(context), ...
  @(x)ischar(x)||isstring(x));
addParameter(p,'IDRSummaryPath',defaultSummaryPath(), ...
  @(x)ischar(x)||isstring(x));
addParameter(p,'TargetPerformance',0.75, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>.5&&x<1);
addParameter(p,'GainBounds',[-5 5], ...
  @(x)isnumeric(x)&&numel(x)==2&&all(isfinite(x))&&x(1)<x(2));
addParameter(p,'SaveOutputs',true,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});opts=p.Results;

dataFolder=analysisPath(context,'Common Code','Data','Simulation');
matPath=fullfile(dataFolder,'IDRIDQ_MTDOGExpectedGains.mat');
csvPath=fullfile(dataFolder,'IDRIDQ_MTDOGExpectedGains.csv');
validateSimulationSourceFiles(context,string(opts.IDRSummaryPath));
if isfile(matPath)
  F=load(matPath,'expected');expected=F.expected;
  if string(expected.context.RunID)~=context.RunID
    error('fitIDRExpectedChoiceDOGGains:RunMismatch', ...
      'Stored expected gains belong to another run.');
  end
  if opts.SaveOutputs&&~isfile(csvPath),writetable(expected.summaryTable,csvPath);end
  fprintf('Loaded completed expected-choice gain fit: %s\n',matPath);
  return
end

if ~isfile(opts.SimulationPath)
  error('fitIDRExpectedChoiceDOGGains:MissingSimulation', ...
    'Synthetic-choice file was not found: %s',opts.SimulationPath);
end
R=load(opts.IDRSummaryPath,'inventory');S=load(opts.SimulationPath,'simulation');
if ~isfield(R,'inventory')||~isfield(S,'simulation')
  error('fitIDRExpectedChoiceDOGGains:InvalidInput', ...
    'Expected inventory and simulation variables.');
end
inventory=R.inventory;simulation=S.simulation;
if string(simulation.context.RunID)~=context.RunID
  error('fitIDRExpectedChoiceDOGGains:SimulationRunMismatch', ...
    'Synthetic choices belong to another run.');
end
if ~isfield(inventory,'sideGainData')|| ...
    ~isfield(inventory,'alignedPsychometric')
  error('fitIDRExpectedChoiceDOGGains:MissingIDRData', ...
    'IDR inventory lacks sideGainData or alignedPsychometric.');
end

D=makeExpectedData(inventory.sideGainData,simulation);
psych=inventory.alignedPsychometric.primary;
syntheticFit=analyzeIDRMatchedGains(D,psych.betaWeibull,psych.lapse, ...
  opts.TargetPerformance,opts.GainBounds);
if isfield(inventory,'matchedGainFits')
  observedFit=inventory.matchedGainFits;
else
  observedFit=analyzeIDRMatchedGains(inventory.sideGainData, ...
    psych.betaWeibull,psych.lapse,opts.TargetPerformance,opts.GainBounds);
end
summaryTable=combineTables(observedFit.summaryTable,syntheticFit.summaryTable);
gCQChiSquare=sum(summaryTable.gCQChiSquareComponent,'omitnan');

expected=struct();expected.createdAt=datetime('now');
expected.createdBy=mfilename;expected.dataOrigin="synthetic_mean_choice_probability";
expected.context=context;expected.nChoiceDraws=simulation.nReplicates;
expected.sigmaReadoutDeg=simulation.sigmaReadoutDeg;
expected.surroundMagnitude=simulation.surroundMagnitude;
expected.sigmaSurroundDeg=simulation.sigmaSurroundDeg;
expected.fitTarget="IDR gCQ(theta) only";
expected.selectionCriterion="sum(((synthetic gCQ-observed gCQ)/observed seCQ)^2)";
expected.gCQChiSquare=gCQChiSquare;expected.summaryTable=summaryTable;
expected.syntheticFit=syntheticFit;

if opts.SaveOutputs
  if ~isfolder(dataFolder),mkdir(dataFolder);end
  expected.outputPaths=string({matPath;csvPath});info=whos('expected');
  assertSimulationRunCapacity(context,info.bytes);
  temporaryMAT=matPath+".temporary.mat";save(temporaryMAT,'expected','-v7.3');
  commitTemporary(temporaryMAT,matPath);
  temporaryCSV=csvPath+".temporary.csv";writetable(summaryTable,temporaryCSV);
  commitTemporary(temporaryCSV,csvPath);assertSimulationRunCapacity(context,0);
end
end

function D=makeExpectedData(source,simulation)
T=source.trialTable;idx=find(cellfun(@(x)string(x.dataset)=="IDR", ...
  simulation.choiceResults));
if isempty(idx)
  error('fitIDRExpectedChoiceDOGGains:NoIDRConditions', ...
    'No IDR choice conditions were found.');
end
session=[];trial=[];offset=[];probability=[];experimental=[];
for k=1:numel(idx)
  C=simulation.choiceResults{idx(k)};
  session=[session;double(C.sessionIndex(:))]; %#ok<AGROW>
  trial=[trial;double(C.trialIndex(:))]; %#ok<AGROW>
  offset=[offset;repmat(double(C.offsetDeg),numel(C.trialIndex),1)]; %#ok<AGROW>
  probability=[probability;double(C.trialProbabilityCorrect(:))]; %#ok<AGROW>
  experimental=[experimental;logical(C.experimentalCorrect(:))]; %#ok<AGROW>
end
row=matchSyntheticTrialRows(T.sessionIndex,T.trialIdx,T.probeDirDeg, ...
  session,trial,offset);
if ~isequal(logical(T.correct(row)),experimental)
  error('fitIDRExpectedChoiceDOGGains:TrialMismatch', ...
    'Mapped IDR trials do not match the simulation inventory.');
end
D=struct();D.trialTable=T(row,:);D.trialTable.correct=probability;
D.stepFrames=source.stepFrames;
fields={'changePrefNoiseByFrameTrial','changeProbeEffectiveNoiseByFrameTrial', ...
  'noChangePrefNoiseByFrameTrial','noChangeProbeEffectiveNoiseByFrameTrial'};
for k=1:numel(fields)
  X=source.(fields{k});D.(fields{k})=X(:,row);
end
end

function T=combineTables(O,S)
[tf,row]=ismember(double(O.probeDirDeg),double(S.probeDirDeg));
if ~all(tf),error('fitIDRExpectedChoiceDOGGains:OffsetMismatch', ...
    'Observed and synthetic gain offsets differ.');end
S=S(row,:);T=table(double(O.probeDirDeg),'VariableNames',{'probeDirDeg'});
names={'gCP','gCQ','gNP','gNQ'};
for k=1:numel(names)
  name=names{k};T.(['observed' name])=O.(name);
  T.(['observedSE' name])=O.(['se' name(2:end)]);
  T.(['synthetic' name])=S.(name);T.(['syntheticSE' name])=S.(['se' name(2:end)]);
end
T.gCQResidual=T.syntheticgCQ-T.observedgCQ;
T.gCQStandardizedResidual=T.gCQResidual./T.observedSEgCQ;
T.gCQChiSquareComponent=T.gCQStandardizedResidual.^2;
end

function commitTemporary(temporary,final)
[ok,message]=movefile(temporary,final,'f');
if ~ok,error('fitIDRExpectedChoiceDOGGains:OutputCommitFailed', ...
    'Could not commit output: %s',message);end
end
function path=defaultSimulationPath(context)
path=analysisPath(context,'Common Code','Data','Simulation', ...
  'IDRIDQ_MTUpstreamNoiseChoices.mat');
end
function path=defaultSummaryPath()
path=fullfile(domainFolder(mfilename('fullpath'),'IDR'), ...
  'Data','AcrossSessionSummaries','IDR_SideGainAnalysis.mat');
end
