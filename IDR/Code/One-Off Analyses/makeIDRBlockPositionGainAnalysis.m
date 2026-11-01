function analysis = makeIDRBlockPositionGainAnalysis(varargin)
% makeIDRBlockPositionGainAnalysis  Test IDR INC gain across block position.
%
% This analysis tests whether probe gain is higher early than late in an INC
% block following a completed DEC block. It consumes the gain-blind block
% inventory produced by makeIDRBlockPositionInventory (version 2 or later),
% and analyzes the fixed-width early, middle, and late windows selected by
% that inventory. Secondary tail trials are not included.
%
% By default, the three requested probe offsets are fit separately. Stable
% differences among session-offsets are removed inside the likelihood by
% fitting one baseline probe gain per parent session. Late is the reference:
%
%   gProbe(trial) = gProbeSession + dProbeEarly*Iearly +
%                                   dProbeMiddle*Imiddle
%
% The primary directional estimand is dProbeEarly = gEarly - gLate.
%
% AnalysisMode='jointLinear' instead fits all offsets together. It assigns
% one baseline probe gain to every observed session-offset and uses a single
% continuous position regressor that equals 1 on the first trial and 0 on
% the final retained trial. The primary estimand is then the common endpoint
% contrast, gProbe(first)-gProbe(last). A secondary model permits a separate
% position slope at each offset as a heterogeneity diagnostic.
%
% Nuisance terms allow block position to alter preferred-noise gain, step
% sensitivity, and side bias. A stable side-bias parameter is also fit for
% each parent session. Parent psychometric alpha, preferred gain, shared beta,
% and shared lapse come from IDRGainAcrossSessionPrefFit_<Animal>.mat.
%
% Confidence intervals use a two-stage cluster bootstrap: parent sessions are
% resampled first, then complete retained INC blocks are resampled intact
% within each selected session. Trials are never resampled individually.
%
% Meetz and Neesha must be run separately. Meetz is the prespecified primary
% animal; Neesha is an independent consistency analysis.
%
% Example
% -------
%   A = makeIDRBlockPositionGainAnalysis('Animal','Meetz');
%   A = makeIDRBlockPositionGainAnalysis('Animal','Neesha');
%   A = makeIDRBlockPositionGainAnalysis('Animal','Meetz', ...
%         'NBootstrap',500,'Replace',true);
%   A = makeIDRBlockPositionGainAnalysis('Animal','Meetz', ...
%         'AnalysisMode','jointLinear','NBootstrap',200,'Replace',true);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'Animal','Meetz',@(x) ischar(x) || isstring(x));
addParameter(p,'AnalysisMode','categorical',@(x) ischar(x) || isstring(x));
addParameter(p,'ProbeDirsDeg',[45 90 135], ...
  @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
addParameter(p,'InventoryPath','',@(x) ischar(x) || isstring(x));
addParameter(p,'NBootstrap',200, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && mod(x,1)==0);
addParameter(p,'BootstrapSeed',314159, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p,'MinTrialsPerBin',100, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(p,'RequireAllBinsPerSession',true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p,'GainLimit',5, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p,'DeltaGainLimit',3, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p,'BiasLimit',3, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p,'StepLogLimit',1.5, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p,'Replace',false,@(x) islogical(x) && isscalar(x));
addParameter(p,'Verbose',true,@(x) islogical(x) && isscalar(x));
parse(p,varargin{:});
opts = p.Results;
opts.ProbeDirsDeg = unique(double(opts.ProbeDirsDeg(:)'));
opts.AnalysisMode = validatestring(opts.AnalysisMode,{'categorical','jointLinear'});

animalTag = char(string(opts.Animal));
if strcmpi(animalTag,'All')
  error('makeIDRBlockPositionGainAnalysis:SeparateAnimals', ...
    'Run Meetz and Neesha separately; Animal=''All'' is not supported.');
end

root = domainFolder(mfilename('fullpath'));
dataRoot = fullfile(root,'Data');
gainFolder = validFolder(fullfile(dataRoot,'AcrossOffsetSummaries','GainAnalysis'));
plotFolder = validFolder(fullfile(root,'Plots','AcrossProbes','GainAnalysis'));

if strlength(string(opts.InventoryPath)) == 0
  inventoryPath = fullfile(gainFolder, ...
    sprintf('IDRBlockPositionInventory_%s.mat',animalTag));
else
  inventoryPath = char(opts.InventoryPath);
end
weightPath = fullfile(dataRoot,'AcrossOffsetSummaries', ...
  sprintf('BetaWeights_%s.mat',animalTag));
prefFitPath = fullfile(gainFolder, ...
  sprintf('IDRGainAcrossSessionPrefFit_%s.mat',animalTag));

if strcmp(opts.AnalysisMode,'jointLinear')
  stem = sprintf('IDRJointLinearBlockPositionGainAnalysis_%s',animalTag);
else
  stem = sprintf('IDRBlockPositionGainAnalysis_%s',animalTag);
end
matPath = fullfile(gainFolder,[stem '.mat']);
if strcmp(opts.AnalysisMode,'jointLinear')
  csvPath = fullfile(gainFolder,[stem '_summary.csv']);
  sessionCSVPath = fullfile(gainFolder,[stem '_cellBaselines.csv']);
else
  csvPath = fullfile(gainFolder,[stem '_offsetSummary.csv']);
  sessionCSVPath = fullfile(gainFolder,[stem '_sessionBaselines.csv']);
end
pdfPath = fullfile(plotFolder,[stem '.pdf']);

if isfile(matPath) && ~opts.Replace
  S = load(matPath,'analysis');
  analysis = S.analysis;
  if opts.Verbose
    fprintf('IDR block-position gain analysis already exists; loaded %s\n',matPath);
  end
  return;
end

requiredFiles = {inventoryPath,weightPath,prefFitPath};
for k = 1:numel(requiredFiles)
  if ~isfile(requiredFiles{k})
    error('makeIDRBlockPositionGainAnalysis:MissingFile', ...
      'Required analysis file not found: %s',requiredFiles{k});
  end
end

I = load(inventoryPath,'inventory');
W = load(weightPath,'weightData');
P = load(prefFitPath,'fitData');
inventory = I.inventory;
weightData = W.weightData;
prefFit = P.fitData;
validateInventory(inventory,animalTag);

if opts.Verbose
  fprintf('Building fixed-window predictors for %s without using tail trials...\n',animalTag);
end
[trialTable,buildSessionTable] = buildTrialTable(inventory,weightData,prefFit,opts);
if isempty(trialTable)
  error('makeIDRBlockPositionGainAnalysis:NoTrials', ...
    'No eligible primary-window INC noise trials were found.');
end

if strcmp(opts.AnalysisMode,'jointLinear')
  analysis = runJointLinearAnalysis(trialTable,buildSessionTable,inventory,opts, ...
    animalTag,inventoryPath,weightPath,prefFitPath,matPath,csvPath, ...
    sessionCSVPath,pdfPath);
  return;
end

offsets = opts.ProbeDirsDeg(:);
rows = repmat(emptyResultRow(),numel(offsets),1);
fits = cell(numel(offsets),1);
bootstrapFits = cell(numel(offsets),1);
sessionBaselineTables = cell(numel(offsets),1);

for iOffset = 1:numel(offsets)
  offset = offsets(iOffset);
  T = trialTable(trialTable.probeDirDeg == offset,:);
  [T,sessionFilterTable] = retainComparableSessions(T,opts);
  if isempty(T)
    row = emptyResultRow();
    row.probeDirDeg = offset;
    row.message = "No comparable sessions.";
    rows(iOffset) = row;
    continue;
  end

  D = tableToData(T);
  counts = accumarray(D.positionBin,1,[3 1]);
  if any(counts < opts.MinTrialsPerBin)
    row = emptyResultRow();
    row.probeDirDeg = offset;
    row.nEarly = counts(1);
    row.nMiddle = counts(2);
    row.nLate = counts(3);
    row.message = sprintf('Too few trials in at least one bin: [%d %d %d].',counts);
    rows(iOffset) = row;
    continue;
  end

  fit = fitPositionModel(D,opts,2);
  B = bootstrapPositionModel(D,opts,offset);
  row = makeResultRow(offset,D,fit,B);
  rows(iOffset) = row;
  fits{iOffset} = fit;
  bootstrapFits{iOffset} = B;
  sessionBaselineTables{iOffset} = makeSessionBaselineTable(T,fit,offset,sessionFilterTable);

  if opts.Verbose
    fprintf(['  %3g deg: n=%d (%d/%d/%d), sessions=%d, blocks=%d, ' ...
      'gE/gM/gL %.3f/%.3f/%.3f, E-L %.3f [%0.3f %0.3f]\n'], ...
      offset,row.nTrials,row.nEarly,row.nMiddle,row.nLate,row.nSessions, ...
      row.nBlocks,row.gEarly,row.gMiddle,row.gLate,row.deltaEarlyLate, ...
      row.deltaEarlyLateCILow,row.deltaEarlyLateCIHigh);
  end
end

resultTable = struct2table(rows);
sessionBaselineTables = sessionBaselineTables(~cellfun(@isempty,sessionBaselineTables));
if isempty(sessionBaselineTables)
  sessionBaselineTable = table();
else
  sessionBaselineTable = vertcat(sessionBaselineTables{:});
end

metadata = struct();
metadata.version = 1;
metadata.analysisName = 'IDRBlockPositionGainAnalysis';
metadata.primaryAnimal = 'Meetz';
metadata.animalRole = conditional(strcmpi(animalTag,'Meetz'), ...
  'primary inference','independent consistency analysis');
metadata.primaryHypothesis = ...
  'Within each offset, normalized probe gain is higher early than late in INC blocks.';
metadata.primaryEstimand = ...
  'dProbeEarly = gProbeEarly - gProbeLate, with one baseline probe gain per session-offset';
metadata.positionWindows = inventory.boundaryTable;
metadata.tailPolicy = 'inventory positionBin==tail excluded from primary analysis';
metadata.sessionNormalization = ...
  'one fitted baseline probe gain per parent session within each offset';
metadata.nuisanceModel = ['early/middle changes in preferred-noise gain, ' ...
  'step sensitivity, and side bias; stable side bias per parent session'];
metadata.psychometricModel = ['fixed parent alpha and preferred gain plus shared beta/lapse; ' ...
  'rectified effective coherence passed through Weibull'];
metadata.bootstrap = sprintf(['%d two-stage resamples: parent sessions, then intact ' ...
  'retained INC blocks within session'],opts.NBootstrap);
metadata.inventoryPath = inventoryPath;
metadata.weightPath = weightPath;
metadata.prefFitPath = prefFitPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

analysis = struct();
analysis.version = 1;
analysis.metadata = metadata;
analysis.inventoryBoundaryTable = inventory.boundaryTable;
analysis.buildSessionTable = buildSessionTable;
analysis.trialTable = trialTable;
analysis.resultTable = resultTable;
analysis.sessionBaselineTable = sessionBaselineTable;
analysis.fits = fits;
analysis.bootstrapFits = bootstrapFits;

save(matPath,'analysis','-v7.3');
writetable(resultTable,csvPath);
if ~isempty(sessionBaselineTable), writetable(sessionBaselineTable,sessionCSVPath); end
makeSummaryFigure(resultTable,pdfPath,animalTag);

if opts.Verbose
  fprintf('Saved MAT: %s\n',matPath);
  fprintf('Saved offset summary: %s\n',csvPath);
  if ~isempty(sessionBaselineTable), fprintf('Saved session baselines: %s\n',sessionCSVPath); end
  fprintf('Saved summary plot: %s\n',pdfPath);
end
end

% =========================================================================
function analysis = runJointLinearAnalysis(T,buildSessionTable,I,opts,animalTag, ...
  inventoryPath,weightPath,prefFitPath,matPath,csvPath,cellCSVPath,pdfPath)

if numel(opts.ProbeDirsDeg) < 2
  error('makeIDRBlockPositionGainAnalysis:JointNeedsOffsets', ...
    'jointLinear mode requires at least two probe offsets.');
end
boundary = I.boundaryTable(strcmpi(string(I.boundaryTable.animal),animalTag),:);
if height(boundary) ~= 1
  error('makeIDRBlockPositionGainAnalysis:JointBoundary', ...
    'Could not identify one inventory boundary row for %s.',animalTag);
end
lastPosition = double(boundary.shortestBlockAllTrials(1));
if lastPosition < 2
  error('makeIDRBlockPositionGainAnalysis:JointBoundary', ...
    'The retained position window must contain at least two trials.');
end
T.positionWeight = (lastPosition-double(T.blockTrialIndex))/(lastPosition-1);
if any(T.positionWeight < -1e-10 | T.positionWeight > 1+1e-10)
  error('makeIDRBlockPositionGainAnalysis:JointPosition', ...
    'A retained trial lies outside the inventory position window.');
end
T.positionWeight = min(max(T.positionWeight,0),1);
[T,cellFilterTable] = retainComparableCells(T,opts);
if isempty(T)
  error('makeIDRBlockPositionGainAnalysis:NoJointCells', ...
    'No session-offset cells contain trials in all three position bins.');
end
D = jointTableToData(T,opts.ProbeDirsDeg);
if numel(unique(D.offsetIndex)) ~= numel(opts.ProbeDirsDeg)
  error('makeIDRBlockPositionGainAnalysis:MissingJointOffset', ...
    'At least one requested offset has no retained session-offset cells.');
end
counts = accumarray(D.positionBin,1,[3 1]);
if any(counts < opts.MinTrialsPerBin)
  error('makeIDRBlockPositionGainAnalysis:JointTooFewTrials', ...
    'Joint position-bin counts are [%d %d %d].',counts);
end

if opts.Verbose
  fprintf(['Fitting joint linear-position model: %d trials, %d sessions, ' ...
    '%d session-offsets, %d blocks...\n'],numel(D.correct), ...
    max(D.sessionIndex),max(D.cellIndex),numel(unique(D.blockKey)));
end
commonFit = fitJointLinearModel(D,opts,'common',2);
offsetFit = fitJointLinearModel(D,opts,'offset',2);
B = bootstrapJointLinear(D,opts);
summaryTable = makeJointSummaryTable(D,commonFit,offsetFit,B,animalTag,opts);
cellBaselineTable = makeJointCellBaselineTable(T,D,commonFit,offsetFit,opts);

metadata = struct();
metadata.version = 1;
metadata.analysisName = 'IDRJointLinearBlockPositionGainAnalysis';
metadata.primaryAnimal = 'Meetz';
metadata.animalRole = conditional(strcmpi(animalTag,'Meetz'), ...
  'primary inference','independent consistency analysis');
metadata.primaryHypothesis = ['A common probe-gain component decreases ' ...
  'continuously from the first to the final retained trial across 45-135 deg offsets.'];
metadata.primaryEstimand = ['dProbePosition = gProbe(first retained trial) - ' ...
  'gProbe(final retained trial), shared across offsets'];
metadata.positionRegressor = sprintf( ...
  'x=(%d-blockTrialIndex)/(%d-1): x=1 first, x=0 final',lastPosition,lastPosition);
metadata.sessionNormalization = 'one fitted probe-gain baseline per session-offset';
metadata.nuisanceModel = ['common linear position terms for preferred-noise gain, ' ...
  'log step sensitivity, and side bias; stable side bias per parent session'];
metadata.secondaryModel = 'separate probe position slope at each offset; diagnostic only';
metadata.tailPolicy = 'inventory tail trials excluded';
metadata.bootstrap = sprintf(['%d two-stage resamples: parent sessions, then intact ' ...
  'retained INC blocks jointly across offsets'],opts.NBootstrap);
metadata.inventoryPath = inventoryPath;
metadata.weightPath = weightPath;
metadata.prefFitPath = prefFitPath;
metadata.options = opts;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

analysis = struct();
analysis.version = 1;
analysis.metadata = metadata;
analysis.inventoryBoundaryTable = boundary;
analysis.buildSessionTable = buildSessionTable;
analysis.cellFilterTable = cellFilterTable;
analysis.trialTable = T;
analysis.summaryTable = summaryTable;
analysis.cellBaselineTable = cellBaselineTable;
analysis.commonFit = commonFit;
analysis.offsetFit = offsetFit;
analysis.bootstrap = B;

save(matPath,'analysis','-v7.3');
writetable(summaryTable,csvPath);
writetable(cellBaselineTable,cellCSVPath);
makeJointSummaryFigure(summaryTable,cellBaselineTable,D,commonFit,offsetFit, ...
  B,pdfPath,animalTag,opts.ProbeDirsDeg);

if opts.Verbose
  fprintf(['  Common early-late slope: %.3f [%0.3f %0.3f], directional p %.4f\n'], ...
    commonFit.decoded.dProbeCommon,B.commonCI(1),B.commonCI(2),B.directionalP);
  fprintf('  Secondary slopes by offset: %s\n', ...
    strjoin(compose('%g: %.3f',opts.ProbeDirsDeg(:), ...
    offsetFit.decoded.dProbeByOffset(:)),', '));
  fprintf('  Heterogeneity DeltaNLL %.3f, approximate p %.4f\n', ...
    commonFit.nll-offsetFit.nll,summaryTable.heterogeneityApproxP);
  fprintf('Saved MAT: %s\n',matPath);
  fprintf('Saved joint summary: %s\n',csvPath);
  fprintf('Saved cell baselines: %s\n',cellCSVPath);
  fprintf('Saved summary plot: %s\n',pdfPath);
end
end

% =========================================================================
function [T,F] = retainComparableCells(T,opts)
pairs = unique([double(T.inventorySessionIndex),double(T.probeDirDeg)],'rows','stable');
rows = cell(size(pairs,1),1);
keep = false(size(pairs,1),1);
for k = 1:size(pairs,1)
  mask = T.inventorySessionIndex==pairs(k,1) & T.probeDirDeg==pairs(k,2);
  counts = accumarray(T.positionBin(mask),1,[3 1]);
  keep(k) = ~opts.RequireAllBinsPerSession || all(counts>0);
  rows{k} = table(pairs(k,1),string(T.sessionName(find(mask,1))),pairs(k,2), ...
    counts(1),counts(2),counts(3),keep(k), ...
    'VariableNames',{'inventorySessionIndex','sessionName','probeDirDeg', ...
    'nEarly','nMiddle','nLate','included'});
end
F = vertcat(rows{:});
keepKey = string(pairs(keep,1)) + "#" + string(pairs(keep,2));
trialKey = string(T.inventorySessionIndex) + "#" + string(T.probeDirDeg);
T = T(ismember(trialKey,keepKey),:);
end

% =========================================================================
function D = jointTableToData(T,offsets)
D = struct();
vars = {'correct','changeSide','deltaC','alpha','gPref','beta','lapse', ...
  'iPref','iProbe','positionBin','positionWeight','blockIndex', ...
  'inventorySessionIndex','probeDirDeg'};
for k = 1:numel(vars), D.(vars{k}) = double(T.(vars{k})(:)); end
sessionValues = unique(D.inventorySessionIndex,'stable');
[tf,D.sessionIndex] = ismember(D.inventorySessionIndex,sessionValues);
if ~all(tf), error('Joint session indexing failed.'); end
D.sessionValues = sessionValues;
D.sessionNames = strings(numel(sessionValues),1);
for k = 1:numel(sessionValues)
  D.sessionNames(k) = T.sessionName(find(D.inventorySessionIndex==sessionValues(k),1));
end
D.offsetValues = double(offsets(:));
D.offsetIndex = zeros(size(D.probeDirDeg));
for k = 1:numel(D.offsetValues)
  D.offsetIndex(abs(D.probeDirDeg-D.offsetValues(k))<1e-9) = k;
end
if any(D.offsetIndex==0), error('Joint offset indexing failed.'); end
cellPairs = unique([D.sessionIndex,D.offsetIndex],'rows','stable');
[tf,D.cellIndex] = ismember([D.sessionIndex,D.offsetIndex],cellPairs,'rows');
if ~all(tf), error('Joint session-offset indexing failed.'); end
D.cellSessionIndex = cellPairs(:,1);
D.cellOffsetIndex = cellPairs(:,2);
D.blockKey = string(T.sessionName) + "#" + string(T.blockIndex);
end

% =========================================================================
function fit = fitJointLinearModel(D,opts,model,nStarts)
if nargin < 4, nStarts=1; end
nCells = max(D.cellIndex);
nSessions = max(D.sessionIndex);
nOffsets = numel(D.offsetValues);
if strcmp(model,'common'), nProbeTerms=1; else, nProbeTerms=nOffsets; end
nParams = nCells+nSessions+nProbeTerms+3;
starts = zeros(nStarts,nParams);
if nStarts>1
  oldRNG=rng; cleanup=onCleanup(@() rng(oldRNG)); %#ok<NASGU>
  rng(424242,'twister'); starts(2:end,:)=0.03*randn(nStarts-1,nParams);
end
objective = @(theta) jointLinearNLL(theta,D,opts,model,nCells,nSessions,nOffsets);
bestNLL=Inf; best=[];
for k=1:nStarts
  if exist('fminunc','file')==2
    o=optimoptions('fminunc','Display','off','Algorithm','quasi-newton', ...
      'MaxIterations',1800,'MaxFunctionEvaluations',150000, ...
      'OptimalityTolerance',1e-7,'StepTolerance',1e-8);
    [theta,nll,exitflag,output]=fminunc(objective,starts(k,:),o);
  else
    o=optimset('Display','off','MaxIter',5000,'MaxFunEvals',200000, ...
      'TolX',1e-7,'TolFun',1e-7);
    [theta,nll,exitflag,output]=fminsearch(objective,starts(k,:),o);
  end
  if isfinite(nll) && nll<bestNLL
    bestNLL=nll; best=struct('theta',theta,'nll',nll, ...
      'exitflag',exitflag,'output',output);
  end
end
if isempty(best), error('All joint optimization starts failed.'); end
best.decoded = decodeJointLinear(best.theta,opts,model,nCells,nSessions,nOffsets);
best.model = model;
fit = best;
end

% =========================================================================
function nll = jointLinearNLL(theta,D,opts,model,nCells,nSessions,nOffsets)
Q = decodeJointLinear(theta,opts,model,nCells,nSessions,nOffsets);
x = D.positionWeight;
gProbe = Q.gProbeCell(D.cellIndex) + Q.dProbeByOffset(D.offsetIndex).*x;
gPref = D.gPref + Q.dPrefPosition.*x;
stepScale = exp(Q.stepLogPosition.*x);
effective = stepScale.*D.deltaC + gPref.*D.iPref + gProbe.*D.iProbe;
effective = max(effective,0);
p0 = 0.5+(0.5-D.lapse).*(1-exp(-((effective./D.alpha).^D.beta)));
p0 = min(max(p0,1e-10),1-1e-10);
sideSign=ones(size(D.changeSide)); sideSign(D.changeSide==1)=-1;
bias = Q.biasSession(D.sessionIndex)+Q.dBiasPosition.*x;
p = stableLogistic(log(p0./(1-p0))+sideSign.*bias);
p = min(max(p,1e-12),1-1e-12);
nll = -sum(D.correct.*log(p)+(1-D.correct).*log(1-p));
if ~isfinite(nll), nll=realmax; end
end

% =========================================================================
function Q = decodeJointLinear(theta,opts,model,nCells,nSessions,nOffsets)
theta=theta(:); k=0;
Q.gProbeCell=opts.GainLimit*tanh(theta(k+(1:nCells))); k=k+nCells;
Q.biasSession=opts.BiasLimit*tanh(theta(k+(1:nSessions))); k=k+nSessions;
if strcmp(model,'common')
  Q.dProbeCommon=opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
  Q.dProbeByOffset=repmat(Q.dProbeCommon,nOffsets,1);
else
  Q.dProbeByOffset=opts.DeltaGainLimit*tanh(theta(k+(1:nOffsets))); k=k+nOffsets;
  Q.dProbeCommon=NaN;
end
Q.dPrefPosition=opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
Q.stepLogPosition=opts.StepLogLimit*tanh(theta(k+1)); k=k+1;
Q.dBiasPosition=opts.BiasLimit*tanh(theta(k+1));
end

% =========================================================================
function B = bootstrapJointLinear(D,opts)
B=struct(); B.nRequested=opts.NBootstrap; B.nSuccessful=0;
B.commonSlope=nan(opts.NBootstrap,1); B.commonCI=[NaN NaN]; B.directionalP=NaN;
if opts.NBootstrap<=0, return; end
oldRNG=rng; cleanup=onCleanup(@() rng(oldRNG)); %#ok<NASGU>
rng(opts.BootstrapSeed+7919,'twister');
for iBoot=1:opts.NBootstrap
  Db=resampleJointSessionsAndBlocks(D);
  try
    fit=fitJointLinearModel(Db,opts,'common',1);
    if isfinite(fit.nll) && isfinite(fit.decoded.dProbeCommon)
      B.commonSlope(iBoot)=fit.decoded.dProbeCommon;
    end
  catch
    % Failed resamples remain NaN and are counted explicitly.
  end
  if opts.Verbose && mod(iBoot,max(1,round(opts.NBootstrap/10)))==0
    fprintf('    joint bootstrap %d/%d\n',iBoot,opts.NBootstrap);
  end
end
ok=isfinite(B.commonSlope); B.nSuccessful=sum(ok);
if B.nSuccessful>=max(20,ceil(0.5*opts.NBootstrap))
  B.commonCI=percentileLocal(B.commonSlope(ok),[2.5 97.5]);
  B.directionalP=(1+sum(B.commonSlope(ok)<=0))/(B.nSuccessful+1);
end
end

% =========================================================================
function Db = resampleJointSessionsAndBlocks(D)
sourceSessions=unique(D.sessionIndex); nSessions=numel(sourceSessions);
sampled=sourceSessions(randi(nSessions,nSessions,1));
indices=cell(0,1); newSession=cell(0,1);
for k=1:nSessions
  mask=D.sessionIndex==sampled(k); blocks=unique(D.blockKey(mask));
  sampledBlocks=blocks(randi(numel(blocks),numel(blocks),1));
  for b=1:numel(sampledBlocks)
    idx=find(mask & D.blockKey==sampledBlocks(b));
    indices{end+1,1}=idx; %#ok<AGROW>
    newSession{end+1,1}=repmat(k,numel(idx),1); %#ok<AGROW>
  end
end
idx=vertcat(indices{:}); Db=subsetData(D,idx);
Db.sessionIndex=vertcat(newSession{:});
cellPairs=unique([Db.sessionIndex,Db.offsetIndex],'rows','stable');
[tf,Db.cellIndex]=ismember([Db.sessionIndex,Db.offsetIndex],cellPairs,'rows');
if ~all(tf), error('Bootstrap cell indexing failed.'); end
Db.cellSessionIndex=cellPairs(:,1); Db.cellOffsetIndex=cellPairs(:,2);
Db.sessionValues=(1:nSessions)'; Db.sessionNames="bootstrap"+string((1:nSessions)');
Db.blockKey="bootstrapTrial"+string((1:numel(idx))');
end

% =========================================================================
function S = makeJointSummaryTable(D,commonFit,offsetFit,B,animalTag,opts)
S=table(); S.animal=string(animalTag); S.nTrials=numel(D.correct);
S.nSessions=max(D.sessionIndex); S.nSessionOffsets=max(D.cellIndex);
S.nBlocks=numel(unique(D.blockKey));
S.commonDeltaEarlyLate=commonFit.decoded.dProbeCommon;
S.commonCILow=B.commonCI(1); S.commonCIHigh=B.commonCI(2);
S.directionalP=B.directionalP; S.bootstrapSuccessful=B.nSuccessful;
S.dPrefEarlyLate=commonFit.decoded.dPrefPosition;
S.stepScaleEarlyVsLate=exp(commonFit.decoded.stepLogPosition);
S.dBiasEarlyLate=commonFit.decoded.dBiasPosition;
S.commonNLL=commonFit.nll; S.commonExitflag=commonFit.exitflag;
S.offsetNLL=offsetFit.nll; S.offsetExitflag=offsetFit.exitflag;
S.heterogeneityDeltaNLL=commonFit.nll-offsetFit.nll;
S.heterogeneityApproxP=gammainc(max(S.heterogeneityDeltaNLL,0), ...
  (numel(opts.ProbeDirsDeg)-1)/2,'upper');
for k=1:numel(opts.ProbeDirsDeg)
  name=matlab.lang.makeValidName(sprintf('secondaryDelta_%gdeg',opts.ProbeDirsDeg(k)));
  S.(name)=offsetFit.decoded.dProbeByOffset(k);
  name=matlab.lang.makeValidName(sprintf('nTrials_%gdeg',opts.ProbeDirsDeg(k)));
  S.(name)=sum(D.offsetIndex==k);
  name=matlab.lang.makeValidName(sprintf('nSessionOffsets_%gdeg',opts.ProbeDirsDeg(k)));
  S.(name)=sum(D.cellOffsetIndex==k);
end
end

% =========================================================================
function C = makeJointCellBaselineTable(T,D,commonFit,offsetFit,opts)
nCells=max(D.cellIndex); C=table();
C.inventorySessionIndex=zeros(nCells,1); C.sessionName=strings(nCells,1);
C.probeDirDeg=zeros(nCells,1); C.nTrials=zeros(nCells,1);
C.gLateCommon=commonFit.decoded.gProbeCell(:);
C.gEarlyCommon=C.gLateCommon+commonFit.decoded.dProbeCommon;
C.gLateOffsetModel=offsetFit.decoded.gProbeCell(:);
C.gEarlyOffsetModel=nan(nCells,1);
for c=1:nCells
  mask=D.cellIndex==c; first=find(mask,1);
  C.inventorySessionIndex(c)=D.inventorySessionIndex(first);
  C.sessionName(c)=T.sessionName(find(T.inventorySessionIndex==D.inventorySessionIndex(first),1));
  C.probeDirDeg(c)=D.probeDirDeg(first); C.nTrials(c)=sum(mask);
  k=find(abs(opts.ProbeDirsDeg-C.probeDirDeg(c))<1e-9,1);
  C.gEarlyOffsetModel(c)=C.gLateOffsetModel(c)+offsetFit.decoded.dProbeByOffset(k);
end
end

% =========================================================================
function makeJointSummaryFigure(S,C,D,commonFit,offsetFit,B,pdfPath,animalTag,offsets)
fig=figure('Color','w','Position',[100 100 1100 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
late=nan(numel(offsets),1);
for k=1:numel(offsets), late(k)=mean(C.gLateCommon(C.probeDirDeg==offsets(k))); end
nexttile; plot(offsets,late+commonFit.decoded.dProbeCommon,'o-','LineWidth',1.3); hold on;
plot(offsets,late,'s-','LineWidth',1.3); yline(0,'k:');
xlabel('Probe offset (deg)'); ylabel('Mean session-offset probe gain');
legend({'First retained trial','Final retained trial'},'Location','best','Box','off'); box off;
nexttile;
errorbar(1,S.commonDeltaEarlyLate,S.commonDeltaEarlyLate-S.commonCILow, ...
  S.commonCIHigh-S.commonDeltaEarlyLate,'o','MarkerFaceColor',[0.1 0.45 0.8], ...
  'Color',[0.1 0.45 0.8],'LineWidth',1.4); hold on; yline(0,'k:');
xlim([0.5 1.5]); xticks(1); xticklabels({'Common'}); ylabel('First - final probe gain'); box off;
nexttile; plot(offsets,offsetFit.decoded.dProbeByOffset,'o-','LineWidth',1.3, ...
  'MarkerFaceColor',[0.85 0.4 0.1]); hold on; yline(0,'k:');
xlabel('Probe offset (deg)'); ylabel('Secondary offset-specific slope'); box off;
nexttile;
counts=zeros(numel(offsets),3);
for k=1:numel(offsets)
  for b=1:3, counts(k,b)=sum(D.offsetIndex==k & D.positionBin==b); end
end
bar(categorical(string(offsets)),counts,'grouped'); xlabel('Probe offset (deg)');
ylabel('Eligible noise trials'); legend({'Early','Middle','Late'},'Box','off'); box off;
sgtitle(sprintf('IDR joint linear block-position gain: %s',animalTag));
exportgraphics(fig,pdfPath,'ContentType','vector'); close(fig);
end

% =========================================================================
function validateInventory(I,animalTag)
required = {'version','metadata','sessionTable','trialTable','boundaryTable'};
missing = required(~isfield(I,required));
if ~isempty(missing)
  error('makeIDRBlockPositionGainAnalysis:BadInventory', ...
    'Inventory is missing: %s',strjoin(missing,', '));
end
trialVars = I.trialTable.Properties.VariableNames;
if I.version < 2 || ~all(ismember({'includePrimaryWindow','blockTrialIndex'},trialVars))
  error('makeIDRBlockPositionGainAnalysis:OldInventory', ...
    'Rebuild the inventory with makeIDRBlockPositionInventory version 2 or later.');
end
animals = unique(string(I.boundaryTable.animal));
if numel(animals) ~= 1 || ~strcmpi(animals,animalTag)
  error('makeIDRBlockPositionGainAnalysis:AnimalMismatch', ...
    'Inventory animal (%s) does not match requested animal (%s).', ...
    strjoin(animals,','),animalTag);
end
end

% =========================================================================
function [T,ST] = buildTrialTable(I,W,PF,opts)

IT = I.trialTable;
sessionRows = repmat(emptyBuildSessionRow(),height(I.sessionTable),1);
trialTables = cell(height(I.sessionTable),1);
nUsed = 0;
weightNames = baseNames(W.sessionFileNames);
fitNames = baseNames(PF.sessionTable.fileName);

for iRow = 1:height(I.sessionTable)
  s = emptyBuildSessionRow();
  s.inventorySessionIndex = I.sessionTable.sessionIndex(iRow);
  s.sessionName = string(I.sessionTable.sessionName(iRow));
  s.filePath = string(I.sessionTable.filePath(iRow));
  if ~strcmpi(string(I.sessionTable.animal(iRow)),string(opts.Animal)) || ...
      I.sessionTable.nRetainedIncBlocks(iRow) < 1
    s.message = "No retained primary blocks for requested animal.";
    sessionRows(iRow) = s;
    continue;
  end

  filePath = char(s.filePath);
  if ~isfile(filePath)
    s.message = "Parent session file not found.";
    sessionRows(iRow) = s;
    continue;
  end
  sessionBase = lower(stripExtension(s.sessionName));
  weightRow = find(weightNames == sessionBase);
  fitRow = find(fitNames == sessionBase);
  s.weightRow = scalarOrNaN(weightRow);
  s.prefFitRow = scalarOrNaN(fitRow);
  if numel(weightRow) ~= 1 || numel(fitRow) ~= 1
    s.message = sprintf('Weight matches=%d, preferred-fit matches=%d.', ...
      numel(weightRow),numel(fitRow));
    sessionRows(iRow) = s;
    continue;
  end

  L = load(filePath,'trials','sessionHeader');
  if ~isfield(L,'trials') || ~isfield(L,'sessionHeader')
    s.message = "Missing trials or sessionHeader.";
    sessionRows(iRow) = s;
    continue;
  end

  invMask = IT.sessionIndex == s.inventorySessionIndex & IT.includePrimaryWindow;
  invRows = IT(invMask,:);
  if isempty(invRows)
    s.message = "No primary-window inventory trials.";
    sessionRows(iRow) = s;
    continue;
  end
  nParent = numel(L.trials);
  useParent = false(nParent,1);
  binByParent = zeros(nParent,1);
  blockByParent = nan(nParent,1);
  blockTrialByParent = nan(nParent,1);
  parentIdx = double(invRows.parentTrialIndex);
  useParent(parentIdx) = true;
  binByParent(parentIdx) = positionBinCode(invRows.positionBin);
  blockByParent(parentIdx) = double(invRows.blockIndex);
  blockTrialByParent(parentIdx) = double(invRows.blockTrialIndex);

  weights = double(W.leaveOneOutWeights(weightRow,:));
  stepTMS = double(W.stepTMS(:)');
  if abs(sum(weights)-1) > 1e-12
    error('makeIDRBlockPositionGainAnalysis:BadWeights', ...
      'Leave-one-out weights do not sum to one for %s.',s.sessionName);
  end

  H = L.sessionHeader;
  frameRateHz = headerScalar(H,'frameRateHz');
  preStepMS = headerScalar(H,'preStepMS');
  stepMS = headerScalar(H,'stepMS');
  frameMS = 1000/frameRateHz;
  nFrames = round((preStepMS+stepMS)/frameMS);
  tMS = (0:nFrames-1)*frameMS;
  stepMask = tMS >= preStepMS & tMS < preStepMS+stepMS;
  if numel(stepTMS) ~= sum(stepMask) || max(abs(stepTMS-tMS(stepMask))) > 1e-6
    error('makeIDRBlockPositionGainAnalysis:WeightTimeMismatch', ...
      'BetaWeights step times do not match %s.',s.sessionName);
  end

  alpha = double(PF.sessionTable.alphaSession(fitRow));
  gPref = double(PF.sessionTable.gPrefSession(fitRow));
  beta = double(PF.params.betaShared);
  lapse = double(PF.params.lapseShared);

  rows = repmat(emptyTrialRow(),sum(useParent),1);
  nRows = 0;
  for iTrial = find(useParent(:)')
    tr = L.trials{iTrial};
    if ~eligibleIncNoiseTrial(tr), continue; end
    D = finalDataRecord(tr.trial.data);
    probeDir = double(D.probeDirDeg);
    if ~any(abs(opts.ProbeDirsDeg-probeDir) < 1e-9), continue; end

    nProbeStreams = 1 + double(probeDir < 180);
    times = numericData(tr.changeTimesMS);
    prefValues = numericData(tr.changePrefCohsPC);
    probeValues = nProbeStreams .* numericData(tr.changeProbeCohsPC);
    [prefFrames,probeFrames] = displayedFrameSequences(times,prefValues, ...
      probeValues,nFrames,frameMS);
    iPref = sum(weights .* prefFrames(stepMask));
    iProbe = sum(weights .* probeFrames(stepMask));

    nRows = nRows+1;
    r = emptyTrialRow();
    r.inventorySessionIndex = s.inventorySessionIndex;
    r.sessionName = s.sessionName;
    r.parentTrialIndex = iTrial;
    r.blockIndex = blockByParent(iTrial);
    r.blockTrialIndex = blockTrialByParent(iTrial);
    r.positionBin = binByParent(iTrial);
    r.probeDirDeg = probeDir;
    r.correct = double(numericData(tr.trialEnd)==0);
    r.changeSide = double(D.changeSide);
    r.deltaC = double(D.stepCohPC)-double(D.baseCohPC);
    r.alpha = alpha;
    r.gPref = gPref;
    r.beta = beta;
    r.lapse = lapse;
    r.iPref = iPref;
    r.iProbe = iProbe;
    rows(nRows) = r;
  end

  rows = rows(1:nRows);
  if nRows > 0
    nUsed = nUsed+1;
    trialTables{nUsed} = struct2table(rows);
    s.included = true;
    s.nTrials = nRows;
    s.message = "ok";
  else
    s.message = "No eligible fixed-window INC noise trials.";
  end
  sessionRows(iRow) = s;
end

trialTables = trialTables(1:nUsed);
if isempty(trialTables), T = table(); else, T = vertcat(trialTables{:}); end
ST = struct2table(sessionRows);
end

% =========================================================================
function [T,S] = retainComparableSessions(T,opts)
if isempty(T)
  S = table();
  return;
end
sessions = unique(T.inventorySessionIndex);
rows = cell(numel(sessions),1);
keepSession = false(numel(sessions),1);
for k = 1:numel(sessions)
  mask = T.inventorySessionIndex == sessions(k);
  counts = accumarray(T.positionBin(mask),1,[3 1]);
  keepSession(k) = ~opts.RequireAllBinsPerSession || all(counts > 0);
  rows{k} = table(sessions(k),string(T.sessionName(find(mask,1))), ...
    counts(1),counts(2),counts(3),keepSession(k), ...
    'VariableNames',{'inventorySessionIndex','sessionName','nEarly', ...
    'nMiddle','nLate','included'});
end
S = vertcat(rows{:});
T = T(ismember(T.inventorySessionIndex,sessions(keepSession)),:);
end

% =========================================================================
function D = tableToData(T)
D = struct();
vars = {'correct','changeSide','deltaC','alpha','gPref','beta','lapse', ...
  'iPref','iProbe','positionBin','blockIndex','inventorySessionIndex'};
for k = 1:numel(vars), D.(vars{k}) = double(T.(vars{k})(:)); end
[sessionValues,~,D.sessionIndex] = unique(D.inventorySessionIndex,'stable');
D.sessionValues = sessionValues;
D.sessionNames = strings(numel(sessionValues),1);
for k = 1:numel(sessionValues)
  D.sessionNames(k) = T.sessionName(find(D.inventorySessionIndex==sessionValues(k),1));
end
D.blockKey = string(T.sessionName) + "#" + string(T.blockIndex);
end

% =========================================================================
function fit = fitPositionModel(D,opts,nStarts)
if nargin < 3, nStarts = 1; end
nSessions = max(D.sessionIndex);
nParams = 2*nSessions + 8;
starts = zeros(nStarts,nParams);
if nStarts > 1
  oldRNG = rng;
  cleanup = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
  rng(8675309,'twister');
  starts(2:end,:) = 0.03*randn(nStarts-1,nParams);
end

objective = @(theta) modelNLL(theta,D,opts,nSessions);
bestNLL = Inf;
best = [];
for k = 1:nStarts
  if exist('fminunc','file') == 2
    o = optimoptions('fminunc','Display','off','Algorithm','quasi-newton', ...
      'MaxIterations',1500,'MaxFunctionEvaluations',100000, ...
      'OptimalityTolerance',1e-7,'StepTolerance',1e-8);
    [theta,nll,exitflag,output] = fminunc(objective,starts(k,:),o);
  else
    o = optimset('Display','off','MaxIter',4000,'MaxFunEvals',150000, ...
      'TolX',1e-7,'TolFun',1e-7);
    [theta,nll,exitflag,output] = fminsearch(objective,starts(k,:),o);
  end
  if isfinite(nll) && nll < bestNLL
    bestNLL = nll;
    best = struct('theta',theta,'nll',nll,'exitflag',exitflag,'output',output);
  end
end
if isempty(best)
  error('makeIDRBlockPositionGainAnalysis:FitFailed','All optimization starts failed.');
end
Q = decodeParameters(best.theta,nSessions,opts);
fit = best;
fit.decoded = Q;
fit.nSessions = nSessions;
fit.meanSessionProbeGain = mean(Q.gProbeSession);
fit.gBins = [fit.meanSessionProbeGain+Q.dProbeEarly, ...
  fit.meanSessionProbeGain+Q.dProbeMiddle,fit.meanSessionProbeGain];
fit.deltaEarlyLate = Q.dProbeEarly;
fit.deltaMiddleLate = Q.dProbeMiddle;
fit.stepScaleBins = [exp(Q.stepLogEarly),exp(Q.stepLogMiddle),1];
fit.prefGainDeltaBins = [Q.dPrefEarly,Q.dPrefMiddle,0];
fit.biasDeltaBins = [Q.dBiasEarly,Q.dBiasMiddle,0];
end

% =========================================================================
function nll = modelNLL(theta,D,opts,nSessions)
Q = decodeParameters(theta,nSessions,opts);
early = D.positionBin==1;
middle = D.positionBin==2;

gProbe = Q.gProbeSession(D.sessionIndex) + ...
  Q.dProbeEarly.*early + Q.dProbeMiddle.*middle;
gPref = D.gPref + Q.dPrefEarly.*early + Q.dPrefMiddle.*middle;
stepScale = exp(Q.stepLogEarly.*early + Q.stepLogMiddle.*middle);
effective = stepScale.*D.deltaC + gPref.*D.iPref + gProbe.*D.iProbe;
effective = max(effective,0);

p0 = 0.5 + (0.5-D.lapse).*(1-exp(-((effective./D.alpha).^D.beta)));
p0 = min(max(p0,1e-10),1-1e-10);

sideSign = ones(size(D.changeSide));
sideSign(D.changeSide==1) = -1;
bias = Q.biasSession(D.sessionIndex) + ...
  Q.dBiasEarly.*early + Q.dBiasMiddle.*middle;
eta = log(p0./(1-p0)) + sideSign.*bias;
p = stableLogistic(eta);
p = min(max(p,1e-12),1-1e-12);
nll = -sum(D.correct.*log(p)+(1-D.correct).*log(1-p));
if ~isfinite(nll), nll = realmax; end
end

% =========================================================================
function Q = decodeParameters(theta,nSessions,opts)
theta = theta(:);
k = 0;
Q.gProbeSession = opts.GainLimit*tanh(theta(k+(1:nSessions))); k=k+nSessions;
Q.biasSession = opts.BiasLimit*tanh(theta(k+(1:nSessions))); k=k+nSessions;
Q.dProbeEarly = opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
Q.dProbeMiddle = opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
Q.dPrefEarly = opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
Q.dPrefMiddle = opts.DeltaGainLimit*tanh(theta(k+1)); k=k+1;
Q.stepLogEarly = opts.StepLogLimit*tanh(theta(k+1)); k=k+1;
Q.stepLogMiddle = opts.StepLogLimit*tanh(theta(k+1)); k=k+1;
Q.dBiasEarly = opts.BiasLimit*tanh(theta(k+1)); k=k+1;
Q.dBiasMiddle = opts.BiasLimit*tanh(theta(k+1));
end

% =========================================================================
function B = bootstrapPositionModel(D,opts,offset)
% Assign fields individually so the cell-valued params field remains inside
% one scalar struct.  Supplying cell(N,1) directly to struct(...) creates an
% N-by-1 struct array (and an empty struct when N is zero).
B = struct();
B.nRequested = opts.NBootstrap;
B.nSuccessful = 0;
B.gBins = nan(opts.NBootstrap,3);
B.deltaEarlyLate = nan(opts.NBootstrap,1);
B.deltaMiddleLate = nan(opts.NBootstrap,1);
B.params = cell(opts.NBootstrap,1);
B.gBinCI = nan(3,2);
B.deltaEarlyLateCI = [NaN NaN];
B.deltaMiddleLateCI = [NaN NaN];
B.directionalP = NaN;
if opts.NBootstrap <= 0, return; end

oldRNG = rng;
cleanup = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
rng(opts.BootstrapSeed+round(offset),'twister');
for iBoot = 1:opts.NBootstrap
  Db = resampleSessionsAndBlocks(D);
  try
    fit = fitPositionModel(Db,opts,1);
    if isfinite(fit.nll) && all(isfinite(fit.gBins))
      B.gBins(iBoot,:) = fit.gBins;
      B.deltaEarlyLate(iBoot) = fit.deltaEarlyLate;
      B.deltaMiddleLate(iBoot) = fit.deltaMiddleLate;
      B.params{iBoot} = fit.decoded;
    end
  catch
    % Leave failed resample as NaN; successful count is reported explicitly.
  end
  if opts.Verbose && mod(iBoot,max(1,round(opts.NBootstrap/10)))==0
    fprintf('    %g deg bootstrap %d/%d\n',offset,iBoot,opts.NBootstrap);
  end
end

ok = isfinite(B.deltaEarlyLate) & all(isfinite(B.gBins),2);
B.nSuccessful = sum(ok);
if B.nSuccessful >= max(20,ceil(0.5*opts.NBootstrap))
  for k = 1:3
    B.gBinCI(k,:) = percentileLocal(B.gBins(ok,k),[2.5 97.5]);
  end
  B.deltaEarlyLateCI = percentileLocal(B.deltaEarlyLate(ok),[2.5 97.5]);
  B.deltaMiddleLateCI = percentileLocal(B.deltaMiddleLate(ok),[2.5 97.5]);
  B.directionalP = (1+sum(B.deltaEarlyLate(ok)<=0))/(B.nSuccessful+1);
end
end

% =========================================================================
function Db = resampleSessionsAndBlocks(D)
sourceSessions = unique(D.sessionIndex);
nSessions = numel(sourceSessions);
sampledSessions = sourceSessions(randi(nSessions,nSessions,1));
indices = cell(0,1);
newSession = cell(0,1);
instance = 0;
for k = 1:numel(sampledSessions)
  source = sampledSessions(k);
  sourceMask = D.sessionIndex==source;
  blocks = unique(D.blockKey(sourceMask));
  sampledBlocks = blocks(randi(numel(blocks),numel(blocks),1));
  instance = instance+1;
  for b = 1:numel(sampledBlocks)
    idx = find(sourceMask & D.blockKey==sampledBlocks(b));
    indices{end+1,1} = idx; %#ok<AGROW>
    newSession{end+1,1} = repmat(instance,numel(idx),1); %#ok<AGROW>
  end
end
idx = vertcat(indices{:});
sessionIndex = vertcat(newSession{:});
Db = subsetData(D,idx);
Db.sessionIndex = sessionIndex;
Db.sessionValues = (1:nSessions)';
Db.sessionNames = "bootstrap"+string((1:nSessions)');
Db.blockKey = "b"+string((1:numel(idx))'); % not used by the fitted likelihood
end

% =========================================================================
function out = subsetData(D,idx)
out = D;
n = numel(D.correct);
fields = fieldnames(D);
for k = 1:numel(fields)
  x = D.(fields{k});
  if (isnumeric(x) || islogical(x) || isstring(x)) && size(x,1)==n
    out.(fields{k}) = x(idx,:);
  end
end
end

% =========================================================================
function row = makeResultRow(offset,D,fit,B)
row = emptyResultRow();
row.probeDirDeg = offset;
row.nTrials = numel(D.correct);
row.nEarly = sum(D.positionBin==1);
row.nMiddle = sum(D.positionBin==2);
row.nLate = sum(D.positionBin==3);
row.nSessions = numel(unique(D.sessionIndex));
row.nBlocks = numel(unique(D.blockKey));
row.correctEarly = mean(D.correct(D.positionBin==1));
row.correctMiddle = mean(D.correct(D.positionBin==2));
row.correctLate = mean(D.correct(D.positionBin==3));
row.gEarly = fit.gBins(1);
row.gMiddle = fit.gBins(2);
row.gLate = fit.gBins(3);
row.gEarlyCILow = B.gBinCI(1,1);
row.gEarlyCIHigh = B.gBinCI(1,2);
row.gMiddleCILow = B.gBinCI(2,1);
row.gMiddleCIHigh = B.gBinCI(2,2);
row.gLateCILow = B.gBinCI(3,1);
row.gLateCIHigh = B.gBinCI(3,2);
row.deltaEarlyLate = fit.deltaEarlyLate;
row.deltaEarlyLateCILow = B.deltaEarlyLateCI(1);
row.deltaEarlyLateCIHigh = B.deltaEarlyLateCI(2);
row.deltaMiddleLate = fit.deltaMiddleLate;
row.deltaMiddleLateCILow = B.deltaMiddleLateCI(1);
row.deltaMiddleLateCIHigh = B.deltaMiddleLateCI(2);
row.directionalP = B.directionalP;
row.dPrefEarly = fit.decoded.dPrefEarly;
row.dPrefMiddle = fit.decoded.dPrefMiddle;
row.stepScaleEarly = fit.stepScaleBins(1);
row.stepScaleMiddle = fit.stepScaleBins(2);
row.dBiasEarly = fit.decoded.dBiasEarly;
row.dBiasMiddle = fit.decoded.dBiasMiddle;
row.nll = fit.nll;
row.exitflag = fit.exitflag;
row.bootstrapSuccessful = B.nSuccessful;
row.message = "ok";
end

% =========================================================================
function S = makeSessionBaselineTable(T,fit,offset,filterTable)
sessions = unique(T.inventorySessionIndex,'stable');
if numel(sessions) ~= numel(fit.decoded.gProbeSession)
  error('makeIDRBlockPositionGainAnalysis:SessionBaselineMismatch', ...
    'Session baseline count does not match fitted session count.');
end
S = table();
S.probeDirDeg = repmat(offset,numel(sessions),1);
S.inventorySessionIndex = sessions;
S.sessionName = strings(numel(sessions),1);
S.gProbeSessionLate = fit.decoded.gProbeSession(:);
S.gProbeSessionEarly = S.gProbeSessionLate + fit.decoded.dProbeEarly;
S.gProbeSessionMiddle = S.gProbeSessionLate + fit.decoded.dProbeMiddle;
S.nEarly = zeros(numel(sessions),1);
S.nMiddle = zeros(numel(sessions),1);
S.nLate = zeros(numel(sessions),1);
for k = 1:numel(sessions)
  mask = T.inventorySessionIndex==sessions(k);
  S.sessionName(k) = T.sessionName(find(mask,1));
  S.nEarly(k) = sum(T.positionBin(mask)==1);
  S.nMiddle(k) = sum(T.positionBin(mask)==2);
  S.nLate(k) = sum(T.positionBin(mask)==3);
end
S.includedAllBins = ismember(S.inventorySessionIndex, ...
  filterTable.inventorySessionIndex(filterTable.included));
end

% =========================================================================
function makeSummaryFigure(T,pdfPath,animalTag)
use = T.message=="ok";
T = T(use,:);
if isempty(T), return; end
fig = figure('Color','w','Position',[100 100 1150 780]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
colors = lines(height(T));

nexttile;
hold on;
for k = 1:height(T)
  x = (1:3) + (k-(height(T)+1)/2)*0.06;
  y = [T.gEarly(k),T.gMiddle(k),T.gLate(k)];
  lo = [T.gEarlyCILow(k),T.gMiddleCILow(k),T.gLateCILow(k)];
  hi = [T.gEarlyCIHigh(k),T.gMiddleCIHigh(k),T.gLateCIHigh(k)];
  errorbar(x,y,y-lo,hi-y,'o-','Color',colors(k,:), ...
    'MarkerFaceColor',colors(k,:),'LineWidth',1.3);
end
yline(0,'k:');
xlim([0.6 3.4]); xticks(1:3); xticklabels({'Early','Middle','Late'});
ylabel('Session-normalized probe gain');
legend(compose('%g deg',T.probeDirDeg),'Location','best','Box','off');
box off;

nexttile;
errorbar(T.probeDirDeg,T.deltaEarlyLate, ...
  T.deltaEarlyLate-T.deltaEarlyLateCILow, ...
  T.deltaEarlyLateCIHigh-T.deltaEarlyLate,'o-', ...
  'Color',[0.15 0.45 0.75],'MarkerFaceColor',[0.15 0.45 0.75], ...
  'LineWidth',1.3);
yline(0,'k:'); xlabel('Probe offset (deg)'); ylabel('Early - late probe gain');
box off;

nexttile;
bar(categorical(string(T.probeDirDeg)),[T.nEarly T.nMiddle T.nLate], ...
  'grouped');
xlabel('Probe offset (deg)'); ylabel('Eligible noise trials');
legend({'Early','Middle','Late'},'Location','best','Box','off'); box off;

nexttile;
hold on;
plot(T.probeDirDeg,T.stepScaleEarly,'o-','LineWidth',1.2);
plot(T.probeDirDeg,T.stepScaleMiddle,'s-','LineWidth',1.2);
yline(1,'k:'); xlabel('Probe offset (deg)'); ylabel('Step-sensitivity scale');
legend({'Early','Middle'},'Location','best','Box','off'); box off;

sgtitle(sprintf('IDR block-position probe gain: %s',animalTag));
exportgraphics(fig,pdfPath,'ContentType','vector');
close(fig);
end

% =========================================================================
function tf = eligibleIncNoiseTrial(tr)
tf = false;
required = {'trial','trialEnd','trialCertify','changeTimesMS', ...
  'changePrefCohsPC','changeProbeCohsPC'};
if ~isstruct(tr) || ~all(isfield(tr,required)) || ...
    ~isstruct(tr.trial) || ~isfield(tr.trial,'data')
  return;
end
D = finalDataRecord(tr.trial.data);
needed = {'changeIndex','cohNoise','probeDirDeg','baseCohPC','stepCohPC', ...
  'changeSide'};
if ~isstruct(D) || ~all(isfield(D,needed)), return; end
outcome = numericData(tr.trialEnd);
certify = numericData(tr.trialCertify);
tf = isscalar(outcome) && ismember(outcome,[0 1]) && ...
  isscalar(certify) && certify==0 && double(D.changeIndex)==1 && ...
  logical(D.cohNoise);
end

% =========================================================================
function [prefFrames,probeFrames] = displayedFrameSequences(times,pref,probe,nFrames,frameMS)
times = double(times(:)); pref = double(pref(:)); probe = double(probe(:));
n = min([numel(times),numel(pref),numel(probe)]);
if n==0
  error('makeIDRBlockPositionGainAnalysis:EmptySequence', ...
    'Eligible noise trial has an empty change-side sequence.');
end
times=times(1:n); pref=pref(1:n); probe=probe(1:n);
if any(diff(times)<0) || times(1)~=0
  error('makeIDRBlockPositionGainAnalysis:BadSequenceTimes', ...
    'Noise event times must begin at zero and be nondecreasing.');
end
eventFrame = floor(times./frameMS)+1;
eventFrame = max(eventFrame,1);
keep = eventFrame<=nFrames;
eventFrame=eventFrame(keep); pref=pref(keep); probe=probe(keep);
prefAtFrame=nan(1,nFrames); probeAtFrame=nan(1,nFrames);
for k=1:numel(eventFrame)
  prefAtFrame(eventFrame(k))=pref(k);
  probeAtFrame(eventFrame(k))=probe(k);
end
eventFrames=find(isfinite(prefAtFrame)&isfinite(probeAtFrame));
if isempty(eventFrames) || eventFrames(1)~=1
  error('makeIDRBlockPositionGainAnalysis:NoInitialState', ...
    'Noise sequence does not define the first displayed frame.');
end
prefFrames=nan(1,nFrames); probeFrames=nan(1,nFrames);
for k=1:numel(eventFrames)
  f0=eventFrames(k);
  if k<numel(eventFrames), f1=eventFrames(k+1)-1; else, f1=nFrames; end
  prefFrames(f0:f1)=prefAtFrame(f0);
  probeFrames(f0:f1)=probeAtFrame(f0);
end
if any(~isfinite(prefFrames)) || any(~isfinite(probeFrames))
  error('makeIDRBlockPositionGainAnalysis:IncompleteSequence', ...
    'Failed to reconstruct every displayed noise frame.');
end
end

% =========================================================================
function code = positionBinCode(labels)
labels = string(labels);
code = zeros(numel(labels),1);
code(labels=="early")=1;
code(labels=="middle")=2;
code(labels=="late")=3;
if any(code==0)
  error('makeIDRBlockPositionGainAnalysis:BadPositionBin', ...
    'Primary-window inventory contains a non-primary position label.');
end
end

% =========================================================================
function x = numericData(x)
while isstruct(x) && isfield(x,'data')
  if numel(x)>1, x=x(end).data; else, x=x.data; end
end
if isnumeric(x) || islogical(x), x=double(x(:)); else, x=[]; end
end

% =========================================================================
function D = finalDataRecord(D)
if isstruct(D) && numel(D)>1, D=D(end); end
end

% =========================================================================
function v = headerScalar(H,name)
if ~isfield(H,name)
  error('makeIDRBlockPositionGainAnalysis:MissingHeader','Missing sessionHeader.%s.',name);
end
v=H.(name);
while isstruct(v) && isfield(v,'data')
  if numel(v)>1, v=v(end).data; else, v=v.data; end
end
v=double(v(1));
end

% =========================================================================
function names = baseNames(x)
x=string(x(:)); names=strings(size(x));
for k=1:numel(x), [~,b,~]=fileparts(x(k)); names(k)=lower(b); end
end

% =========================================================================
function s = stripExtension(x)
[~,b,~]=fileparts(char(string(x))); s=string(b);
end

% =========================================================================
function x = scalarOrNaN(x)
if isscalar(x), x=double(x); else, x=NaN; end
end

% =========================================================================
function p = stableLogistic(x)
p=zeros(size(x)); pos=x>=0;
p(pos)=1./(1+exp(-x(pos)));
z=exp(x(~pos)); p(~pos)=z./(1+z);
end

% =========================================================================
function q = percentileLocal(x,pct)
x=sort(double(x(isfinite(x))));
q=nan(size(pct));
if isempty(x), return; end
for k=1:numel(pct)
  pos=1+(numel(x)-1)*pct(k)/100;
  lo=floor(pos); hi=ceil(pos);
  if lo==hi, q(k)=x(lo); else, q(k)=x(lo)+(pos-lo)*(x(hi)-x(lo)); end
end
end

% =========================================================================
function out = conditional(tf,a,b)
if tf, out=a; else, out=b; end
end

% =========================================================================
function r = emptyBuildSessionRow()
r=struct('inventorySessionIndex',NaN,'sessionName',"",'filePath',"", ...
  'weightRow',NaN,'prefFitRow',NaN,'included',false,'nTrials',0,'message',"");
end

% =========================================================================
function r = emptyTrialRow()
r=struct('inventorySessionIndex',NaN,'sessionName',"", ...
  'parentTrialIndex',NaN,'blockIndex',NaN,'blockTrialIndex',NaN,'positionBin',NaN, ...
  'probeDirDeg',NaN,'correct',NaN,'changeSide',NaN,'deltaC',NaN, ...
  'alpha',NaN,'gPref',NaN,'beta',NaN,'lapse',NaN,'iPref',NaN,'iProbe',NaN);
end

% =========================================================================
function r = emptyResultRow()
r=struct('probeDirDeg',NaN,'nTrials',0,'nEarly',0,'nMiddle',0,'nLate',0, ...
  'nSessions',0,'nBlocks',0,'correctEarly',NaN,'correctMiddle',NaN, ...
  'correctLate',NaN,'gEarly',NaN,'gMiddle',NaN,'gLate',NaN, ...
  'gEarlyCILow',NaN,'gEarlyCIHigh',NaN,'gMiddleCILow',NaN, ...
  'gMiddleCIHigh',NaN,'gLateCILow',NaN,'gLateCIHigh',NaN, ...
  'deltaEarlyLate',NaN,'deltaEarlyLateCILow',NaN, ...
  'deltaEarlyLateCIHigh',NaN,'deltaMiddleLate',NaN, ...
  'deltaMiddleLateCILow',NaN,'deltaMiddleLateCIHigh',NaN, ...
  'directionalP',NaN,'dPrefEarly',NaN,'dPrefMiddle',NaN, ...
  'stepScaleEarly',NaN,'stepScaleMiddle',NaN,'dBiasEarly',NaN, ...
  'dBiasMiddle',NaN,'nll',NaN,'exitflag',NaN,'bootstrapSuccessful',0, ...
  'message',"");
end
