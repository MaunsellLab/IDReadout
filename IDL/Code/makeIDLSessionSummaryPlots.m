function pdfPaths = makeIDLSessionSummaryPlots(varargin)
% makeIDLSessionSummaryPlots  Create one PDF summary per IDL FullSession.
%
% FullSessions are the authoritative inputs selected by selectAnalysisFiles.
% A compact Data/SessionAnalysis file is created when missing and reused on
% later calls. PDFs are saved in IDL/Plots/FullSessions.

p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('TwoAllocationMode',[],@(x) isempty(x)||(islogical(x)&&isscalar(x)));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('ReplaceAnalysis', false, @(x) islogical(x) && isscalar(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;
names=fieldnames(p.Unmatched); selection=cell(1,2*numel(names));
for k=1:numel(names)
  selection{2*k-1}=names{k}; selection{2*k}=p.Unmatched.(names{k});
end
makeIDLSessionAnalyses('Animal',opts.Animal,'TwoAllocationMode',opts.TwoAllocationMode, ...
  'Replace',opts.ReplaceAnalysis,'Verbose',opts.Verbose,selection{:});

root = domainFolder(mfilename('fullpath'));
sourceFolder = fullfile(root, 'Data', 'FullSessions');
sourceFiles = selectAnalysisFiles(sourceFolder,'Animal',opts.Animal, ...
  'TwoAllocationMode',opts.TwoAllocationMode,selection{:});
analysisFolder = validFolder(fullfile(root, 'Data', 'SessionAnalysis'));
plotFolder = idlPlotFolder('FullSessions');

pdfPaths = strings(numel(sourceFiles), 1);
for k = 1:numel(sourceFiles)
  sourcePath = sourceFiles{k};
  [~, sessionBase] = fileparts(sourcePath);
  analysisPath = fullfile(analysisFolder, ...
    [sessionBase '_sessionAnalysis.mat']);
  pdfPath = fullfile(plotFolder, [sessionBase '_IDLSummary.pdf']);
  pdfPaths(k) = string(pdfPath);

  if isfile(analysisPath) && ~opts.ReplaceAnalysis
    loaded = load(analysisPath, 'sessionAnalysis');
    sessionAnalysis = loaded.sessionAnalysis;
  else
    sessionAnalysis = makeIDLSessionAnalysis(sourcePath);
    sessionHeader = sessionAnalysis.sessionHeader; %#ok<NASGU>
    save(analysisPath, 'sessionAnalysis', 'sessionHeader', '-v7.3');
    if opts.Verbose
      fprintf('makeIDLSessionSummaryPlots: saved %s\n', analysisPath);
    end
  end

  if isfile(pdfPath) && ~opts.Replace && ~opts.ReplaceAnalysis
    pdfInfo=dir(pdfPath); analysisInfo=dir(analysisPath);
    if pdfInfo.datenum>=analysisInfo.datenum, continue, end
  end
  [summary, fig] = plotIDLSessionSummary(sessionAnalysis, 'PDFPath', pdfPath, 'Visible', 'off');
  close(fig);
  assert(strcmp(summary.pdfPath, string(pdfPath)), 'makeIDLSessionSummaryPlots:UnexpectedPDFPath', ...
    'The session plot was not written to its expected destination.');
  if opts.Verbose
    fprintf('makeIDLSessionSummaryPlots: saved %s\n', pdfPath);
  end
end
end
