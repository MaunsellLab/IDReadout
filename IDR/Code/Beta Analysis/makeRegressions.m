function batch = makeRegressions(replace, varargin)
% makeRegressions  Produce authoritative per-probe-session beta regressions.
%
% Inputs are Data/Probe*/ProbeSessions/*.mat. Outputs remain in the flat
% Data/Probe*/Regression folder using the historical filename suffix:
%   *_scalarNoiseRegression.mat
%
% Old pilot products are replaced automatically because they do not contain
% the current versioned kernelWeightedProbeRegression structure.

if nargin < 1 || isempty(replace), replace = false; end
p = inputParser;
p.FunctionName = mfilename;
addParameter(p,'MakePlots',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Verbose', false ,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
dataRoot = fullfile(root,'Data');
plotRoot = fullfile(root,'Plots','Probes');
weightPath = fullfile(dataRoot, 'AcrossOffsetSummaries', 'BetaWeights.mat');
W = load(weightPath, 'weightData');
weightData = W.weightData;

probeDirs = dir(fullfile(dataRoot,'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);
batch = struct('version',1,'replace',replace,'weightPath', weightPath,'files',{{}},'regPaths',{{}},'plotPaths',{{}}, ...
  'ok',false(0,1), 'skippedIneligible',false(0,1), 'messages',{{}},'createdBy',mfilename,'createdDate',datetime('now'));
for iProbe = 1:numel(probeDirs)
  probeTag = probeDirs(iProbe).name;
  sessionFolder = fullfile(probeDirs(iProbe).folder,probeTag,'ProbeSessions');
  if ~isfolder(sessionFolder), continue; end
  regressionFolder = validFolder(fullfile(probeDirs(iProbe).folder,probeTag,'Regression'));
  plotFolder = fullfile(plotRoot,probeTag,'Regression');
  if opts.MakePlots, validFolder(plotFolder); end
  files = dir(fullfile(sessionFolder,'*.mat'));
  [~,order] = sort({files.name}); 
  files = files(order);
  for iFile = 1:numel(files)
    sourcePath = fullfile(files(iFile).folder,files(iFile).name);
    [~,baseName] = fileparts(sourcePath);
    regPath = fullfile(regressionFolder,[baseName '_scalarNoiseRegression.mat']);
    plotPath = fullfile(plotFolder,[baseName '_scalarNoiseRegression.pdf']);
    batch.files{end+1,1}=sourcePath; 
    batch.regPaths{end+1,1}=regPath; 
    batch.plotPaths{end+1,1}=plotPath; 
    batch.ok(end+1,1)=false; 
    batch.skippedIneligible(end+1,1)=false; 
    batch.messages{end+1,1}=''; 
    ib = numel(batch.ok);

    current = isCurrentProduct(regPath);
    needReg = replace || ~current;
    needPlot = opts.MakePlots && (replace || needReg || ~isfile(plotPath));
    if ~needReg && ~needPlot
      if opts.Verbose, fprintf('Skipping current regression: %s [%s]\n',baseName,probeTag); end
      batch.ok(ib)=true; batch.messages{ib}='current outputs; skipped';
      continue;
    end

    % Skip probe sessions whose parent full session was not included when BetaWeights.mat was constructed.
    Shead = load(sourcePath, 'sessionHeader', 'sessionProbeHeader');
    parentName = '';
    if isfield(Shead, 'sessionProbeHeader') && ...
        isfield(Shead.sessionProbeHeader, 'parentFileName')
      parentName = char(Shead.sessionProbeHeader.parentFileName);
    elseif isfield(Shead, 'sessionHeader') && ...
        isfield(Shead.sessionHeader, 'fileName')
      parentName = char(Shead.sessionHeader.fileName);
    end

    [~, parentBase] = fileparts(parentName);
    weightBases = cellfun(@(x) filepartsBase(x), ...
      weightData.sessionFileNames, 'UniformOutput', false);

    nWeightMatches = sum(strcmpi(parentBase, weightBases));

    if nWeightMatches == 0
      if opts.Verbose
        fprintf('Skipping ineligible session: %s [%s] — parent absent from BetaWeights\n', ...
          baseName, probeTag);
      end
      batch.messages{ib} = 'parent session absent from BetaWeights; skipped';
      batch.skippedIneligible(ib) = true;
      continue;
    elseif nWeightMatches > 1
      warning('makeRegressions:DuplicateWeightMatch', ...
        'Skipping %s: parent session matched %d BetaWeights rows.', ...
        sourcePath, nWeightMatches);
      batch.messages{ib} = 'parent session matched multiple BetaWeights rows';
      continue;
    end

    try
      if needReg
        if opts.Verbose, fprintf('Fitting regression: %s [%s]\n',baseName,probeTag); end
        reg = computeKernelWeightedProbeRegression(sourcePath,weightData);
        reg.weightInfo.weightSourceFile = weightPath;
        save(regPath,'reg','-v7.3');
      else
        R = load(regPath,'reg'); reg = R.reg;
      end
      if needPlot
        fig = plotKernelWeightedProbeRegression(reg);
        exportgraphics(fig,plotPath,'ContentType','vector');
        close(fig);
      end
      batch.ok(ib)=true; batch.messages{ib}='ok';
      if opts.Verbose && reg.fitByStep.inc.fitUsable
        fprintf('  INC beta pref %.4g, probe %.4g, ratio %.4g\n', ...
          reg.fitByStep.inc.betaPref,reg.fitByStep.inc.betaProbe,reg.fitByStep.inc.betaRatio);
      end
    catch ME
      warning('makeRegressions:SessionFailed','Failed on %s:\n%s',sourcePath,ME.message);
      batch.messages{ib}=ME.message;
    end
  end
end

summaryPath = fullfile(dataRoot,'scalarNoiseRegression_batchSummary.mat');
save(summaryPath,'batch');
if opts.Verbose
  fprintf('\nRegression production complete: %d/%d successful.\n',sum(batch.ok),numel(batch.ok));
  nProduced = sum(batch.ok);
  nIneligible = sum(batch.skippedIneligible);
  nFailed = numel(batch.ok) - nProduced - nIneligible;
  fprintf('\nRegression production complete:\n');
  fprintf('  successful:          %d\n', nProduced);
  fprintf('  skipped ineligible:  %d\n', nIneligible);
  fprintf('  failed:              %d\n', nFailed);
  fprintf('Batch summary saved: %s\n', summaryPath);
end
end

function tf = isCurrentProduct(path)
tf = false;
if ~isfile(path), return; end
try
  S = load(path,'reg');
  tf = isfield(S,'reg') && isfield(S.reg,'version') && S.reg.version >= 3 && ...
    isfield(S.reg,'analysisName') && strcmp(S.reg.analysisName,'kernelWeightedProbeRegression');
catch
  tf = false;
end
end

function base = filepartsBase(pathOrName)
[~, base] = fileparts(char(pathOrName));
end