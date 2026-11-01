function pdfPaths = makeIDJSessionSummaryPlots(varargin)
% makeIDJSessionSummaryPlots  Create one PDF summary per IDJ FullSession.
%
% FullSessions are the authoritative inputs selected by selectAnalysisFiles.
% A compact Data/SessionAnalysis file is created when missing and reused on
% later calls. PDFs are saved in IDJ/Plots/FullSessions.

p = inputParser;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('ReplaceAnalysis', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

root = domainFolder(mfilename('fullpath'));
sourceFolder = fullfile(root, 'Data', 'FullSessions');
sourceFiles = selectAnalysisFiles(sourceFolder, 'Animal', opts.Animal);
analysisFolder = validFolder(fullfile(root, 'Data', 'SessionAnalysis'));
plotFolder = idjPlotFolder('FullSessions');

pdfPaths = strings(numel(sourceFiles), 1);
for k = 1:numel(sourceFiles)
  sourcePath = sourceFiles{k};
  [~, sessionBase] = fileparts(sourcePath);
  analysisPath = fullfile(analysisFolder, ...
    [sessionBase '_sessionAnalysis.mat']);
  pdfPath = fullfile(plotFolder, [sessionBase '_IDJSummary.pdf']);
  pdfPaths(k) = string(pdfPath);

  analysisCurrent = false;
  if isfile(analysisPath) && ~opts.ReplaceAnalysis
    loaded = load(analysisPath, 'sessionAnalysis');
    if isfield(loaded,'sessionAnalysis') && idjHasTopUpExclusion(loaded.sessionAnalysis)
      sessionAnalysis = loaded.sessionAnalysis;
      analysisCurrent = true;
    end
  end
  if ~analysisCurrent
    sessionAnalysis = makeIDJSessionAnalysis(sourcePath);
    sessionHeader = sessionAnalysis.sessionHeader; %#ok<NASGU>
    save(analysisPath, 'sessionAnalysis', 'sessionHeader', '-v7.3');
    if opts.Verbose
      fprintf('makeIDJSessionSummaryPlots: saved %s\n', analysisPath);
    end
  end

  pdfCurrent = false;
  if isfile(pdfPath)
    pdfInfo = dir(pdfPath);
    analysisInfo = dir(analysisPath);
    pdfCurrent = pdfInfo.datenum >= analysisInfo.datenum;
  end
  if pdfCurrent && analysisCurrent && ~opts.Replace && ~opts.ReplaceAnalysis
    continue
  end
  [summary, fig] = plotIDJSessionSummary(sessionAnalysis, 'PDFPath', pdfPath, 'Visible', 'off');
  close(fig);
  assert(strcmp(summary.pdfPath, string(pdfPath)), 'makeIDJSessionSummaryPlots:UnexpectedPDFPath', ...
    'The session plot was not written to its expected destination.');
  if opts.Verbose
    fprintf('makeIDJSessionSummaryPlots: saved %s\n', pdfPath);
  end
end
end
