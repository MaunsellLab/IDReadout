function [sessionAnalysisAll, rebuilt] = makeIDJSessionAnalyses(varargin)
% makeIDJSessionAnalyses  Build compact analyses for converted IDJ sessions.

p = inputParser;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
sourceFolder = fullfile(root, 'Data', 'FullSessions');
outputFolder = validFolder(fullfile(root, 'Data', 'SessionAnalysis'));
sourceFiles = selectAnalysisFiles(sourceFolder, 'Animal', opts.Animal);

sessionAnalysisAll = cell(numel(sourceFiles), 1);
rebuilt = false(numel(sourceFiles),1);
for k = 1:numel(sourceFiles)
  sourcePath = sourceFiles{k};
  [~, sourceBase] = fileparts(sourcePath);
  outputPath = fullfile(outputFolder, [sourceBase '_sessionAnalysis.mat']);
  if isfile(outputPath) && ~opts.Replace
    vars = who('-file', outputPath);
    if all(ismember({'sessionAnalysis','sessionHeader'}, vars))
      existing = load(outputPath, 'sessionAnalysis');
      if idjHasTopUpExclusion(existing.sessionAnalysis)
        sessionAnalysisAll{k} = existing.sessionAnalysis;
        continue
      end
    end
  end

  sessionAnalysis = makeIDJSessionAnalysis(sourcePath);
  rebuilt(k) = true;
  sessionHeader = sessionAnalysis.sessionHeader; %#ok<NASGU>
  save(outputPath, 'sessionAnalysis', 'sessionHeader', '-v7.3');
  sessionAnalysisAll{k} = sessionAnalysis;
  if opts.Verbose
    fprintf('makeIDJSessionAnalyses: saved %s\n', outputPath);
  end
end

if isscalar(sessionAnalysisAll)
  sessionAnalysisAll = sessionAnalysisAll{1};
end
end
