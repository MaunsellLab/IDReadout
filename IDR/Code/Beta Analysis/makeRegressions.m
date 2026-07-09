function makeRegressions(varargin)
% makeRegressions  Produce authoritative per-probe-session beta regressions.
%
% Inputs are Data/Probe*/ProbeSessions/*.mat. Outputs remain in the flat
% Data/Probe*/Regression folder using the historical filename suffix:
%   *_scalarNoiseRegression.mat
%
% Old pilot products are replaced automatically because they do not contain
% the current versioned kernelWeightedProbeRegression structure.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) isempty(x) || ischar(x) || isstring(x));
addParameter(p, 'MakePlots', true, @(x)islogical(x)&&isscalar(x));
addParameter(p, 'Replace', false , @(x)islogical(x)&&isscalar(x));
addParameter(p, 'Verbose', false , @(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opts = p.Results;

% We use a weighted kernel that should be updated before running
makeBetaKernel('Animal', opts.Animal);

root = domainFolder(mfilename('fullpath'));
dataRoot = fullfile(root, 'Data');
plotRoot = fullfile(root, 'Plots', 'Probes');
weightPath = fullfile(dataRoot, 'AcrossOffsetSummaries', sprintf('BetaWeights_%s.mat', opts.Animal));
W = load(weightPath, 'weightData');
weightData = W.weightData;

probeDirs = dir(fullfile(dataRoot,'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);
for iProbe = 1:numel(probeDirs)
  probeTag = probeDirs(iProbe).name;
  sessionFolder = fullfile(probeDirs(iProbe).folder,probeTag,'ProbeSessions');
  if ~isfolder(sessionFolder), continue; end
  regressionFolder = validFolder(fullfile(probeDirs(iProbe).folder, probeTag, 'Regression'));
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
    current = isCurrentProduct(regPath);
    needReg = p.Results.Replace || ~current;
    needPlot = opts.MakePlots && (p.Results.Replace || needReg || ~isfile(plotPath));
    if ~needReg && ~needPlot
      if opts.Verbose, fprintf('Skipping current regression: %s [%s]\n',baseName,probeTag); end
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
      continue;
    elseif nWeightMatches > 1
      warning('makeRegressions:DuplicateWeightMatch', 'Skipping %s: parent session matched %d BetaWeights rows.', ...
        sourcePath, nWeightMatches);
      continue;
    end

    try
      if needReg
        if opts.Verbose, fprintf('Fitting regression: %s [%s]\n',baseName,probeTag); end
        reg = computeKernelWeightedProbeRegression(sourcePath, weightData);
        reg.weightInfo.weightSourceFile = weightPath;
        sessionHeader = Shead.sessionHeader;
        save(regPath,'reg', 'sessionHeader', '-v7.3');
      else
        R = load(regPath, 'reg'); reg = R.reg;
      end
      if needPlot
        % fig = plotKernelWeightedProbeRegression(reg);
        % exportgraphics(fig, plotPath, 'ContentType', 'vector');
        % close(fig);
      end
      if opts.Verbose && reg.fitByStep.inc.fitUsable
        fprintf('  INC beta pref %.4g, probe %.4g, ratio %.4g\n', ...
          reg.fitByStep.inc.betaPref,reg.fitByStep.inc.betaProbe,reg.fitByStep.inc.betaRatio);
      end
    catch ME
      warning('makeRegressions:SessionFailed','Failed on %s:\n%s',sourcePath,ME.message);
    end
  end
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

