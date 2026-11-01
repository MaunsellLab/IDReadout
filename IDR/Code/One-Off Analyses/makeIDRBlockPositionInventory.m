function inventory = makeIDRBlockPositionInventory(varargin)
% makeIDRBlockPositionInventory  Survey IDR block lengths and eligibility.
%
% This is the gain-blind first stage of the IDR block-position analysis. It
% reconstructs INC/DEC blocks from every parent-session trial, certifies a
% block against the eligible-trial quota stored in blockStatus.repsPerBlock,
% applies the prespecified structural exclusions, and selects fixed
% early/middle/late boundaries from the shortest retained block.
%
% No stimulus-noise predictor or behavioral gain is calculated or loaded.
%
% Primary retained blocks
% -----------------------
% * INC blocks only.
% * Standard noise amplitudes (default 10% preferred, 7.071% probe).
% * At least one requested probe offset (default 45, 90, or 135 deg).
% * First physical block of each session excluded.
% * Incomplete blocks excluded.
% * INC block must immediately follow a complete DEC block.
% * Final N complete INC blocks of each session excluded (default N=2).
%
% Block position counts every parent-session trial, whether eligible or not.
% Eligibility uses the established IDR criterion:
%   trialCertify == 0 and trialEnd is 0 or 1.
%
% Boundaries are chosen separately for each animal. If Lmin is the shortest
% retained complete INC block in all-trial units:
%   early  = 1 : floor(Lmin/3)
%   middle = floor(Lmin/3)+1 : floor(2*Lmin/3)
%   late   = floor(2*Lmin/3)+1 : Lmin
% Trials after Lmin are labeled tail and are not part of the primary thirds.
%
% Example
% -------
%   I = makeIDRBlockPositionInventory('Animal','Meetz');
%   I = makeIDRBlockPositionInventory('Animal','Neesha');
%   I = makeIDRBlockPositionInventory('Animal','All','Replace',true);

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(p, 'TargetProbeDirsDeg', [45 90 135], ...
  @(x) isnumeric(x) && isvector(x) && all(isfinite(x)));
addParameter(p, 'PrefCohNoisePC', 10, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'ProbeCohNoisePC', 5*sqrt(2), ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x > 0);
addParameter(p, 'NoiseTolerancePC', 0.01, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
addParameter(p, 'ExcludeFirstBlock', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'NLastCompleteIncBlocksToExclude', 2, ...
  @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && mod(x,1) == 0);
addParameter(p, 'RequireCompletePrecedingDec', true, ...
  @(x) islogical(x) && isscalar(x));
addParameter(p, 'DataFolder', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'OutputFolder', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'Replace', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;
opts.TargetProbeDirsDeg = unique(double(opts.TargetProbeDirsDeg(:)'));

root = domainFolder(mfilename('fullpath'));
if strlength(string(opts.DataFolder)) == 0
  dataFolder = fullfile(root, 'Data', 'FullSessions');
else
  dataFolder = char(opts.DataFolder);
end
if strlength(string(opts.OutputFolder)) == 0
  outputFolder = validFolder(fullfile(root, 'Data', ...
    'AcrossOffsetSummaries', 'GainAnalysis'));
else
  outputFolder = validFolder(char(opts.OutputFolder));
end

animalTag = char(string(opts.Animal));
stem = sprintf('IDRBlockPositionInventory_%s', animalTag);
matPath = fullfile(outputFolder, [stem '.mat']);
blockCSVPath = fullfile(outputFolder, [stem '_blocks.csv']);
positionCSVPath = fullfile(outputFolder, [stem '_eligibilityByPosition.csv']);
sessionEndCSVPath = fullfile(outputFolder, [stem '_eligibilityBySessionEnd.csv']);
boundaryCSVPath = fullfile(outputFolder, [stem '_boundaries.csv']);
pdfPath = fullfile(outputFolder, [stem '.pdf']);

if isfile(matPath) && ~opts.Replace
  S = load(matPath, 'inventory');
  inventory = S.inventory;
  if opts.Verbose
    fprintf('IDR block-position inventory already exists; loaded %s\n', matPath);
  end
  return;
end

[selectedFiles, selectedInfo] = selectAnalysisFiles(dataFolder, ...
  'Animal', opts.Animal, 'Bin180Into179', false);
if isempty(selectedFiles)
  error('makeIDRBlockPositionInventory:NoSessions', ...
    'No IDR parent-session files passed selectAnalysisFiles.');
end

if opts.Verbose
  fprintf('Surveying %d parent sessions without calculating gains...\n', ...
    numel(selectedFiles));
end

blockTables = cell(numel(selectedFiles),1);
trialTables = cell(numel(selectedFiles),1);
sessionRows = repmat(emptySessionRow(), numel(selectedFiles), 1);

for iFile = 1:numel(selectedFiles)
  filePath = selectedFiles{iFile};
  [~, sessionName] = fileparts(filePath);
  L = load(filePath, 'trials', 'sessionHeader');

  s = emptySessionRow();
  s.sessionIndex = iFile;
  s.sessionName = string(sessionName);
  s.filePath = string(filePath);
  if ~isfield(L, 'trials') || ~isfield(L, 'sessionHeader')
    s.message = "Missing trials or sessionHeader.";
    sessionRows(iFile) = s;
    continue;
  end

  animal = headerText(L.sessionHeader, {'animal','subject'}, '');
  [BT, TT] = inventoryOneSession(L.trials, animal, sessionName, iFile, opts);
  blockTables{iFile} = BT;
  trialTables{iFile} = TT;

  s.animal = string(animal);
  s.nTrials = numel(L.trials);
  s.nBlocks = height(BT);
  s.nIncBlocks = sum(BT.isInc);
  s.nCompleteIncBlocks = sum(BT.isInc & BT.isComplete);
  s.nRetainedIncBlocks = sum(BT.includePrimary);
  s.message = "ok";
  sessionRows(iFile) = s;
end

blockTables = blockTables(~cellfun(@isempty,blockTables));
trialTables = trialTables(~cellfun(@isempty,trialTables));
if isempty(blockTables)
  error('makeIDRBlockPositionInventory:NoUsableSessions', ...
    'No selected file supplied usable trial and block metadata.');
end
blockTable = vertcat(blockTables{:});
trialTable = vertcat(trialTables{:});
sessionTable = struct2table(sessionRows);

% Choose boundaries independently for each animal, using only retained
% all-trial block lengths. Then label every trial in those blocks.
[boundaryTable, trialTable] = assignBoundaries(trialTable, blockTable);
positionTable = summarizeEligibilityByPosition(trialTable);
binTable = summarizeEligibilityByBin(trialTable);
sessionEndTable = summarizeEligibilityBySessionEnd(blockTable);

metadata = struct();
metadata.version = 2;
metadata.analysisName = 'IDRBlockPositionInventory';
metadata.gainBlind = true;
metadata.method = ['Reconstruct contiguous INC/DEC blocks from all parent ' ...
  'trials; certify completion from repsPerBlock and eligible-trial counts; ' ...
  'select fixed thirds from the shortest retained all-trial block per animal.'];
metadata.eligibilityDefinition = ...
  'trialCertify==0 and trialEnd in {0,1}';
metadata.completionDefinition = ...
  'eligible count reaches blockStatus.repsPerBlock quota for block change type';
metadata.boundaryDefinition = ...
  ['earlyEnd=floor(Lmin/3); middleEnd=floor(2*Lmin/3); ' ...
  'late ends at Lmin; trials after Lmin are secondary tail trials'];
metadata.options = opts;
metadata.dataFolder = dataFolder;
metadata.selectedFileInfo = selectedInfo;
metadata.createdBy = mfilename;
metadata.createdDate = datetime('now');

inventory = struct();
inventory.version = 2;
inventory.metadata = metadata;
inventory.sessionTable = sessionTable;
inventory.blockTable = blockTable;
inventory.trialTable = trialTable;
inventory.boundaryTable = boundaryTable;
inventory.positionTable = positionTable;
inventory.binTable = binTable;
inventory.sessionEndTable = sessionEndTable;

save(matPath, 'inventory', '-v7.3');
writetable(blockTable, blockCSVPath);
writetable(positionTable, positionCSVPath);
writetable(sessionEndTable, sessionEndCSVPath);
writetable(boundaryTable, boundaryCSVPath);
makeInventoryFigure(blockTable, positionTable, sessionEndTable, ...
  boundaryTable, pdfPath);

if opts.Verbose
  fprintf('\nIDR block-position inventory (no gains examined):\n');
  disp(boundaryTable);
  fprintf('Saved MAT: %s\n', matPath);
  fprintf('Saved block table: %s\n', blockCSVPath);
  fprintf('Saved position table: %s\n', positionCSVPath);
  fprintf('Saved session-end table: %s\n', sessionEndCSVPath);
  fprintf('Saved boundary table: %s\n', boundaryCSVPath);
  fprintf('Saved summary plot: %s\n', pdfPath);
end
end

% =========================================================================
function [BT, TT] = inventoryOneSession(trials, animal, sessionName, sessionIndex, opts)

nTrials = numel(trials);
changeIndex = nan(nTrials,1);
isEligible = false(nTrials,1);
hasStepNoise = false(nTrials,1);
probeDirDeg = nan(nTrials,1);
prefCohNoisePC = nan(nTrials,1);
probeCohNoisePC = nan(nTrials,1);
requiredEligible = nan(nTrials,1);

for iTrial = 1:nTrials
  tr = trials{iTrial};
  if ~isstruct(tr)
    continue;
  end

  outcome = fieldNumeric(tr, 'trialEnd');
  certify = fieldNumeric(tr, 'trialCertify');
  isEligible(iTrial) = isscalar(outcome) && ismember(outcome,[0 1]) && ...
    isscalar(certify) && certify == 0;

  if isfield(tr,'trial') && isstruct(tr.trial) && isfield(tr.trial,'data')
    D = tr.trial.data;
    % Some older converted event fields can contain more than one record
    % within a trial. The final record is the state in force for that trial.
    if isstruct(D) && numel(D) > 1, D = D(end); end
    if isfield(D,'changeIndex')
      changeIndex(iTrial) = double(D.changeIndex) + 1; % 1=DEC, 2=INC
    end
    if isfield(D,'cohNoise'), hasStepNoise(iTrial) = logical(D.cohNoise); end
    if isfield(D,'probeDirDeg'), probeDirDeg(iTrial) = double(D.probeDirDeg); end
    if isfield(D,'prefCohNoisePC'), prefCohNoisePC(iTrial) = double(D.prefCohNoisePC); end
    if isfield(D,'probeCohNoisePC'), probeCohNoisePC(iTrial) = double(D.probeCohNoisePC); end
  end

  if isfinite(changeIndex(iTrial)) && isfield(tr,'blockStatus') && ...
      isstruct(tr.blockStatus) && isfield(tr.blockStatus,'data')
    B = tr.blockStatus.data;
    % Older sessions sometimes store multiple blockStatus events in one
    % trial. Select the final status instead of allowing B.repsPerBlock to
    % expand as a comma-separated list of function arguments.
    if isstruct(B) && numel(B) > 1, B = B(end); end
    if isstruct(B) && isscalar(B) && isfield(B,'repsPerBlock')
      q = numericData(B.repsPerBlock);
      k = changeIndex(iTrial);
      if numel(q) >= k
        requiredEligible(iTrial) = double(q(k));
      end
    end
  end
end

% Missing task-state labels are not expected, but isolated omissions should
% not split a physical block. Fill from the previous state, then backward
% from the next state for any omissions at the start of the file.
changeIndex = fillTaskState(changeIndex);
if any(~ismember(changeIndex,[1 2]))
  error('makeIDRBlockPositionInventory:BadChangeIndex', ...
    'Could not reconstruct INC/DEC state for every trial in %s.', sessionName);
end

blockStart = [1; find(diff(changeIndex) ~= 0) + 1];
blockEnd = [blockStart(2:end)-1; nTrials];
nBlocks = numel(blockStart);

blockRows = repmat(emptyBlockRow(),nBlocks,1);
blockIndexByTrial = nan(nTrials,1);
blockTrialIndex = nan(nTrials,1);

for iBlock = 1:nBlocks
  idx = (blockStart(iBlock):blockEnd(iBlock))';
  blockIndexByTrial(idx) = iBlock;
  blockTrialIndex(idx) = (1:numel(idx))';
  state = changeIndex(idx(1));

  quotas = requiredEligible(idx);
  quotas = quotas(isfinite(quotas) & quotas >= 0);
  uniqueQuotas = unique(quotas);
  if isempty(uniqueQuotas)
    quota = NaN;
  else
    quota = uniqueQuotas(end);
  end

  noiseIdx = idx(hasStepNoise(idx));
  offsets = unique(probeDirDeg(noiseIdx));
  offsets = offsets(isfinite(offsets) & offsets > 0);
  prefAmps = uniqueRounded(prefCohNoisePC(noiseIdx));
  probeAmps = uniqueRounded(probeCohNoisePC(noiseIdx));

  r = emptyBlockRow();
  r.sessionIndex = sessionIndex;
  r.sessionName = string(sessionName);
  r.animal = string(animal);
  r.blockIndex = iBlock;
  r.changeIndex = state;
  r.changeType = string(conditional(state == 2,'INC','DEC'));
  r.isInc = state == 2;
  r.startTrialIdx = idx(1);
  r.endTrialIdx = idx(end);
  r.nAllTrials = numel(idx);
  r.nEligibleTrials = sum(isEligible(idx));
  r.nIneligibleTrials = r.nAllTrials - r.nEligibleTrials;
  r.eligibilityRate = r.nEligibleTrials / r.nAllTrials;
  r.requiredEligibleTrials = quota;
  r.quotaChangedWithinBlock = numel(uniqueQuotas) > 1;
  r.isComplete = isfinite(quota) && r.nEligibleTrials >= quota;
  r.isFirstPhysicalBlock = iBlock == 1;
  r.isTerminalPhysicalBlock = iBlock == nBlocks;
  r.blocksFromSessionEnd = nBlocks - iBlock;
  r.probeOffsetsDeg = joinNumeric(offsets);
  r.prefCohNoisePC = scalarIfUnique(prefAmps);
  r.probeCohNoisePC = scalarIfUnique(probeAmps);
  r.hasTargetProbeOffset = any(ismembertol(offsets, ...
    opts.TargetProbeDirsDeg, 1e-9, 'DataScale', 1));
  r.hasStandardNoise = isscalar(prefAmps) && isscalar(probeAmps) && ...
    abs(prefAmps - opts.PrefCohNoisePC) <= opts.NoiseTolerancePC && ...
    abs(probeAmps - opts.ProbeCohNoisePC) <= opts.NoiseTolerancePC;
  if iBlock > 1
    r.precededByDec = changeIndex(blockStart(iBlock)-1) == 1;
  end
  blockRows(iBlock) = r;
end

% A preceding DEC block is complete if its own eligible count reaches its
% quota. Mark the last complete INC block before applying exclusions.
for iBlock = 2:nBlocks
  blockRows(iBlock).precededByCompleteDec = ...
    blockRows(iBlock-1).changeIndex == 1 && blockRows(iBlock-1).isComplete;
end
completeInc = find([blockRows.isInc] & [blockRows.isComplete]);
if ~isempty(completeInc)
  blockRows(completeInc(end)).isLastCompleteIncBlock = true;
  for k = 1:numel(completeInc)
    blockRows(completeInc(k)).completeIncBlocksFromEnd = ...
      numel(completeInc) - k;
  end
  nTerminalExclude = min(opts.NLastCompleteIncBlocksToExclude,numel(completeInc));
  if nTerminalExclude > 0
    terminalExclude = completeInc(end-nTerminalExclude+1:end);
    for k = 1:numel(terminalExclude)
      blockRows(terminalExclude(k)).isTerminalExcludedIncBlock = true;
    end
  end
end

for iBlock = 1:nBlocks
  r = blockRows(iBlock);
  reasons = strings(0,1);
  if ~r.isInc, reasons(end+1) = "not INC"; end %#ok<AGROW>
  if opts.ExcludeFirstBlock && r.isFirstPhysicalBlock
    reasons(end+1) = "first physical block"; %#ok<AGROW>
  end
  if ~r.isComplete, reasons(end+1) = "incomplete quota"; end %#ok<AGROW>
  if ~r.hasTargetProbeOffset, reasons(end+1) = "no target probe offset"; end %#ok<AGROW>
  if ~r.hasStandardNoise, reasons(end+1) = "nonstandard noise"; end %#ok<AGROW>
  if opts.RequireCompletePrecedingDec && ~r.precededByCompleteDec
    reasons(end+1) = "not preceded by complete DEC"; %#ok<AGROW>
  end
  if r.isTerminalExcludedIncBlock
    reasons(end+1) = sprintf('final %d complete INC block(s)', ...
      opts.NLastCompleteIncBlocksToExclude); %#ok<AGROW>
  end
  r.includePrimary = isempty(reasons);
  r.exclusionReason = strjoin(reasons,'; ');
  blockRows(iBlock) = r;
end

BT = struct2table(blockRows);

TT = table();
TT.sessionIndex = repmat(sessionIndex,nTrials,1);
TT.sessionName = repmat(string(sessionName),nTrials,1);
TT.animal = repmat(string(animal),nTrials,1);
TT.parentTrialIndex = (1:nTrials)';
TT.blockIndex = blockIndexByTrial;
TT.blockTrialIndex = blockTrialIndex;
TT.changeIndex = changeIndex;
TT.isEligible = isEligible;
TT.hasStepNoise = hasStepNoise;
TT.probeDirDeg = probeDirDeg;
TT.prefCohNoisePC = prefCohNoisePC;
TT.probeCohNoisePC = probeCohNoisePC;
TT.includePrimaryBlock = false(nTrials,1);
for iBlock = 1:nBlocks
  TT.includePrimaryBlock(blockIndexByTrial == iBlock) = ...
    blockRows(iBlock).includePrimary;
end
TT.positionBin = repmat("excluded",nTrials,1);
TT.includePrimaryWindow = false(nTrials,1);
end

% =========================================================================
function [B, T] = assignBoundaries(T, blockTable)

animals = unique(blockTable.animal(blockTable.includePrimary));
rows = repmat(emptyBoundaryRow(),numel(animals),1);

for iAnimal = 1:numel(animals)
  animal = animals(iAnimal);
  useBlocks = blockTable.includePrimary & blockTable.animal == animal;
  lengths = blockTable.nAllTrials(useBlocks);
  if isempty(lengths)
    continue;
  end
  Lmin = min(lengths);
  earlyEnd = floor(Lmin/3);
  middleEnd = floor(2*Lmin/3);
  if earlyEnd < 1 || middleEnd <= earlyEnd
    error('makeIDRBlockPositionInventory:BlocksTooShort', ...
      'Retained blocks for %s are too short to define thirds.', animal);
  end

  useTrials = T.includePrimaryBlock & T.animal == animal;
  early = useTrials & T.blockTrialIndex <= earlyEnd;
  middle = useTrials & T.blockTrialIndex > earlyEnd & ...
    T.blockTrialIndex <= middleEnd;
  late = useTrials & T.blockTrialIndex > middleEnd & ...
    T.blockTrialIndex <= Lmin;
  tail = useTrials & T.blockTrialIndex > Lmin;
  T.positionBin(early) = "early";
  T.positionBin(middle) = "middle";
  T.positionBin(late) = "late";
  T.positionBin(tail) = "tail";
  T.includePrimaryWindow(early | middle | late) = true;

  r = emptyBoundaryRow();
  r.animal = animal;
  r.nRetainedBlocks = sum(useBlocks);
  r.shortestBlockAllTrials = Lmin;
  r.longestBlockAllTrials = max(lengths);
  r.earlyEndAllTrial = earlyEnd;
  r.middleEndAllTrial = middleEnd;
  r.nEarlyAllTrials = sum(early);
  r.nMiddleAllTrials = sum(middle);
  r.nLateAllTrials = sum(late);
  r.nTailAllTrials = sum(tail);
  r.nEarlyEligible = sum(T.isEligible(early));
  r.nMiddleEligible = sum(T.isEligible(middle));
  r.nLateEligible = sum(T.isEligible(late));
  r.nTailEligible = sum(T.isEligible(tail));
  r.earlyEligibilityRate = mean(T.isEligible(early));
  r.middleEligibilityRate = mean(T.isEligible(middle));
  r.lateEligibilityRate = mean(T.isEligible(late));
  r.tailEligibilityRate = mean(T.isEligible(tail));
  rows(iAnimal) = r;
end

B = struct2table(rows);
end

% =========================================================================
function P = summarizeEligibilityByPosition(T)

use = T.includePrimaryBlock;
animals = unique(T.animal(use));
rows = cell(0,1);
for iAnimal = 1:numel(animals)
  animal = animals(iAnimal);
  positions = unique(T.blockTrialIndex(use & T.animal == animal));
  for iPosition = 1:numel(positions)
    pos = positions(iPosition);
    mask = use & T.animal == animal & T.blockTrialIndex == pos;
    r = table(animal,pos,sum(mask),sum(T.isEligible(mask)),mean(T.isEligible(mask)), ...
      'VariableNames',{'animal','blockTrialIndex','nBlocksContributing', ...
      'nEligible','eligibilityRate'});
    rows{end+1,1} = r; %#ok<AGROW>
  end
end
if isempty(rows), P = table(); else, P = vertcat(rows{:}); end
end

% =========================================================================
function B = summarizeEligibilityByBin(T)

use = T.includePrimaryBlock;
animals = unique(T.animal(use));
bins = ["early","middle","late"];
rows = cell(0,1);
for iAnimal = 1:numel(animals)
  for iBin = 1:numel(bins)
    mask = use & T.animal == animals(iAnimal) & T.positionBin == bins(iBin);
    r = table(animals(iAnimal),bins(iBin),sum(mask),sum(T.isEligible(mask)), ...
      mean(T.isEligible(mask)), ...
      'VariableNames',{'animal','positionBin','nAllTrials','nEligible', ...
      'eligibilityRate'});
    rows{end+1,1} = r; %#ok<AGROW>
  end
end
if isempty(rows), B = table(); else, B = vertcat(rows{:}); end
end

% =========================================================================
function E = summarizeEligibilityBySessionEnd(BT)

% This diagnostic deliberately restores the last complete INC block. It
% uses behavior only and is intended to show whether eligibility was already
% deteriorating in the block removed by the conservative primary analysis.
use = BT.isInc & BT.isComplete & ~BT.isFirstPhysicalBlock & ...
  BT.precededByCompleteDec & BT.hasTargetProbeOffset & BT.hasStandardNoise;
animals = unique(BT.animal(use));
rows = cell(0,1);
for iAnimal = 1:numel(animals)
  distances = unique(BT.completeIncBlocksFromEnd(use & BT.animal == animals(iAnimal)));
  distances = sort(distances(:));
  for iDistance = 1:numel(distances)
    distance = distances(iDistance);
    mask = use & BT.animal == animals(iAnimal) & ...
      BT.completeIncBlocksFromEnd == distance;
    r = table(animals(iAnimal),distance,sum(mask), ...
      mean(BT.nAllTrials(mask)),sum(BT.nEligibleTrials(mask)), ...
      sum(BT.nAllTrials(mask)), ...
      sum(BT.nEligibleTrials(mask))/sum(BT.nAllTrials(mask)), ...
      'VariableNames',{'animal','completeIncBlocksFromEnd','nBlocks', ...
      'meanBlockAllTrials','nEligible','nAllTrials','eligibilityRate'});
    rows{end+1,1} = r; %#ok<AGROW>
  end
end
if isempty(rows), E = table(); else, E = vertcat(rows{:}); end
end

% =========================================================================
function makeInventoryFigure(BT, PT, sessionEndTable, boundaryTable, pdfPath)

if isempty(boundaryTable), return; end
if isfile(pdfPath), delete(pdfPath); end

for iAnimal = 1:height(boundaryTable)
  animal = boundaryTable.animal(iAnimal);
  useBlocks = BT.includePrimary & BT.animal == animal;
  P = PT(PT.animal == animal,:);
  End = sessionEndTable(sessionEndTable.animal == animal,:);

  fig = figure('Color','w','Position',[100 100 1100 760]);
  tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

  nexttile;
  histogram(BT.nAllTrials(useBlocks),'BinMethod','integers', ...
    'FaceColor',[0.25 0.5 0.8]);
  hold on;
  xline(boundaryTable.earlyEndAllTrial(iAnimal),'k--','Early end');
  xline(boundaryTable.middleEndAllTrial(iAnimal),'k--','Middle end');
  xline(boundaryTable.shortestBlockAllTrials(iAnimal),'k--','Late end');
  xlabel('All trials per retained INC block');
  ylabel('Number of blocks');
  box off;

  nexttile;
  yyaxis left;
  plot(P.blockTrialIndex,P.eligibilityRate,'k-','LineWidth',1.2);
  hold on;
  xline(boundaryTable.earlyEndAllTrial(iAnimal),'--','Color',[0.3 0.3 0.3]);
  xline(boundaryTable.middleEndAllTrial(iAnimal),'--','Color',[0.3 0.3 0.3]);
  xline(boundaryTable.shortestBlockAllTrials(iAnimal),'--','Color',[0.3 0.3 0.3]);
  ylim([0 1]);
  ylabel('Eligibility rate');
  yyaxis right;
  plot(P.blockTrialIndex,P.nBlocksContributing,'-','Color',[0.65 0.65 0.65], ...
    'LineWidth',0.9);
  ylabel('Blocks contributing');
  xlabel('All-trial position within INC block');
  box off;

  nexttile;
  plot(End.completeIncBlocksFromEnd,End.eligibilityRate,'o-', ...
    'Color',[0.75 0.25 0.2],'MarkerFaceColor',[0.75 0.25 0.2], ...
    'LineWidth',1.2);
  ylim([0 1]);
  ylabel('Eligibility rate');
  xlabel('Complete INC blocks from end (0 = last)');
  box off;

  nexttile;
  axis off;
  summaryText = sprintf([ ...
    '%s\n\nRetained blocks: %d\nShortest/longest: %d / %d all trials\n' ...
    'Early: 1-%d\nMiddle: %d-%d\nLate: %d-%d\nTail: %d-block end\n\n' ...
    'Eligible trials: %d / %d / %d\nEligibility: %.3f / %.3f / %.3f\n' ...
    'Secondary tail: %d eligible, rate %.3f'], ...
    animal, boundaryTable.nRetainedBlocks(iAnimal), ...
    boundaryTable.shortestBlockAllTrials(iAnimal), ...
    boundaryTable.longestBlockAllTrials(iAnimal), ...
    boundaryTable.earlyEndAllTrial(iAnimal), ...
    boundaryTable.earlyEndAllTrial(iAnimal)+1, ...
    boundaryTable.middleEndAllTrial(iAnimal), ...
    boundaryTable.middleEndAllTrial(iAnimal)+1, ...
    boundaryTable.shortestBlockAllTrials(iAnimal), ...
    boundaryTable.shortestBlockAllTrials(iAnimal)+1, ...
    boundaryTable.nEarlyEligible(iAnimal), ...
    boundaryTable.nMiddleEligible(iAnimal), ...
    boundaryTable.nLateEligible(iAnimal), ...
    boundaryTable.earlyEligibilityRate(iAnimal), ...
    boundaryTable.middleEligibilityRate(iAnimal), ...
    boundaryTable.lateEligibilityRate(iAnimal), ...
    boundaryTable.nTailEligible(iAnimal), ...
    boundaryTable.tailEligibilityRate(iAnimal));
  text(0,1,summaryText,'VerticalAlignment','top','FontName','Courier', ...
    'FontSize',11,'Interpreter','none');
  sgtitle(sprintf('IDR block-position inventory: %s (gain blind)',animal));

  exportgraphics(fig,pdfPath,'ContentType','vector','Append',iAnimal > 1);
  close(fig);
end
end

% =========================================================================
function x = fillTaskState(x)
for i = 2:numel(x)
  if ~isfinite(x(i)), x(i) = x(i-1); end
end
for i = numel(x)-1:-1:1
  if ~isfinite(x(i)), x(i) = x(i+1); end
end
end

% =========================================================================
function x = fieldNumeric(S, fieldName)
x = [];
if ~isfield(S,fieldName), return; end
x = numericData(S.(fieldName));
end

% =========================================================================
function x = numericData(x)
while isstruct(x) && isfield(x,'data')
  if numel(x) > 1
    x = x(end).data;
  else
    x = x.data;
  end
end
if isnumeric(x) || islogical(x)
  x = double(x(:));
else
  x = [];
end
end

% =========================================================================
function x = uniqueRounded(x)
x = double(x(:));
x = x(isfinite(x));
x = unique(round(x,6));
end

% =========================================================================
function x = scalarIfUnique(x)
if isscalar(x), x = double(x); else, x = NaN; end
end

% =========================================================================
function s = joinNumeric(x)
if isempty(x)
  s = "";
else
  s = strjoin(compose('%g',x(:)'),',');
end
end

% =========================================================================
function out = conditional(tf,a,b)
if tf, out = a; else, out = b; end
end

% =========================================================================
function value = headerText(H, names, defaultValue)
value = defaultValue;
for i = 1:numel(names)
  if ~isfield(H,names{i}), continue; end
  value = H.(names{i});
  while isstruct(value) && isfield(value,'data'), value = value.data; end
  value = char(string(value));
  return;
end
end

% =========================================================================
function r = emptySessionRow()
r = struct('sessionIndex',NaN,'sessionName',"",'filePath',"",'animal',"", ...
  'nTrials',0,'nBlocks',0,'nIncBlocks',0,'nCompleteIncBlocks',0, ...
  'nRetainedIncBlocks',0,'message',"");
end

% =========================================================================
function r = emptyBlockRow()
r = struct('sessionIndex',NaN,'sessionName',"",'animal',"", ...
  'blockIndex',NaN,'changeIndex',NaN,'changeType',"",'isInc',false, ...
  'startTrialIdx',NaN,'endTrialIdx',NaN,'nAllTrials',0, ...
  'nEligibleTrials',0,'nIneligibleTrials',0,'eligibilityRate',NaN, ...
  'requiredEligibleTrials',NaN,'quotaChangedWithinBlock',false, ...
  'isComplete',false,'isFirstPhysicalBlock',false, ...
  'isTerminalPhysicalBlock',false,'blocksFromSessionEnd',NaN, ...
  'completeIncBlocksFromEnd',NaN, ...
  'probeOffsetsDeg',"",'prefCohNoisePC',NaN,'probeCohNoisePC',NaN, ...
  'hasTargetProbeOffset',false,'hasStandardNoise',false, ...
  'precededByDec',false,'precededByCompleteDec',false, ...
  'isLastCompleteIncBlock',false,'isTerminalExcludedIncBlock',false, ...
  'includePrimary',false, ...
  'exclusionReason',"");
end

% =========================================================================
function r = emptyBoundaryRow()
r = struct('animal',"",'nRetainedBlocks',0,'shortestBlockAllTrials',NaN, ...
  'longestBlockAllTrials',NaN,'earlyEndAllTrial',NaN, ...
  'middleEndAllTrial',NaN,'nEarlyAllTrials',0,'nMiddleAllTrials',0, ...
  'nLateAllTrials',0,'nTailAllTrials',0,'nEarlyEligible',0, ...
  'nMiddleEligible',0,'nLateEligible',0,'nTailEligible',0, ...
  'earlyEligibilityRate',NaN,'middleEligibilityRate',NaN, ...
  'lateEligibilityRate',NaN,'tailEligibilityRate',NaN);
end
