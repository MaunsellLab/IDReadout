function sessionAnalysisAll = makeIDLSessionAnalyses(varargin)
% makeIDLSessionAnalyses  Build compact analyses for converted IDL sessions.

p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('TwoAllocationMode',[],@(x) isempty(x)||(islogical(x)&&isscalar(x)));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
sourceFolder = fullfile(root, 'Data', 'FullSessions');
outputFolder = validFolder(fullfile(root, 'Data', 'SessionAnalysis'));
selection = namedPairs(p.Unmatched);
sourceFiles = selectAnalysisFiles(sourceFolder, 'Animal', opts.Animal,selection{:});
% Upgrade FullSession headers before applying the new fast header filter.
for k=1:numel(sourceFiles), idlRefreshSessionHeader(sourceFiles{k}); end
sourceFiles = selectAnalysisFiles(sourceFolder,'Animal',opts.Animal, ...
  'TwoAllocationMode',opts.TwoAllocationMode,selection{:});

sessionAnalysisAll = cell(numel(sourceFiles), 1);
for k = 1:numel(sourceFiles)
  sourcePath = sourceFiles{k};
  [~, sourceBase] = fileparts(sourcePath);
  outputPath = fullfile(outputFolder, [sourceBase '_sessionAnalysis.mat']);
  if isfile(outputPath) && ~opts.Replace
    vars = who('-file', outputPath);
    if all(ismember({'sessionAnalysis','sessionHeader'}, vars))
      existing = load(outputPath, 'sessionAnalysis');
      A = existing.sessionAnalysis;
      if isfield(A,'version') && A.version>=2 && ...
          isfield(A.sessionHeader,'idlAllocationVersion') && ...
          A.sessionHeader.idlAllocationVersion==1
        sessionAnalysisAll{k} = A;
        continue
      end
    end
  end

  sessionAnalysis = makeIDLSessionAnalysis(sourcePath);
  sessionHeader = sessionAnalysis.sessionHeader; %#ok<NASGU>
  save(outputPath, 'sessionAnalysis', 'sessionHeader', '-v7.3');
  sessionAnalysisAll{k} = sessionAnalysis;
  if opts.Verbose
    fprintf('makeIDLSessionAnalyses: saved %s\n', outputPath);
  end
end


if isscalar(sessionAnalysisAll)
  sessionAnalysisAll = sessionAnalysisAll{1};
end
end

function pairs=namedPairs(S)
names=fieldnames(S); pairs=cell(1,2*numel(names));
for k=1:numel(names), pairs{2*k-1}=names{k}; pairs{2*k}=S.(names{k}); end
end
