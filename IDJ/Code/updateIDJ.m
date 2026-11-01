function updateIDJ(varargin)
% updateIDJ  Convert IDJ sessions, build session summaries and across-session analyses.
% AcrossSessions defaults to true; set false for session summaries only.
% PNormProfile is opt-in because its full-session bootstrap is expensive.
p = inputParser;
p.addParameter('Animal', 'All', @(x) ischar(x) || isstring(x));
p.addParameter('Replace', false, @(x) islogical(x) && isscalar(x));
p.addParameter('ReplaceAnalysis', false, @(x) islogical(x) && isscalar(x));
p.addParameter('AcrossSessions', true, @(x) islogical(x) && isscalar(x));
p.addParameter('PNormProfile', true, @(x) islogical(x) && isscalar(x));
p.addParameter('NBootstrap', 1000, @(x) isnumeric(x) && isscalar(x) && ...
  isfinite(x) && x>=0 && x==round(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;

fprintf('>>> updateIDJ start\n');
convertIDJData();
[~, rebuilt] = makeIDJSessionAnalyses('Animal', opts.Animal, 'Replace', opts.ReplaceAnalysis, 'Verbose', opts.Verbose);
makeIDJSessionSummaryPlots('Animal', opts.Animal, ...
  'Replace', opts.Replace || opts.ReplaceAnalysis || any(rebuilt), 'ReplaceAnalysis', false, 'Verbose', opts.Verbose);
if opts.AcrossSessions
  result = fitIDJModelRecoveryAcrossSessions('Animal', opts.Animal, 'MakePlots', false);
  plotIDJAcrossSessionBehavior('Animal', opts.Animal, 'ModelResult', result);
end
if opts.PNormProfile
  runPNormProfiles(opts);
end
% 
fprintf('<<< updateIDJ complete\n');
end

function runPNormProfiles(opts)
root = domainFolder(mfilename('fullpath'));
paths = selectAnalysisFiles(fullfile(root,'Data','SessionAnalysis'), 'Animal',opts.Animal);
if ischar(paths), paths = {paths}; end
if isstring(paths), paths = cellstr(paths(:)); end
animals = strings(numel(paths),1);
for k = 1:numel(paths)
  H = load(paths{k},'sessionHeader');
  animals(k) = string(H.sessionHeader.animal);
end
for animal = unique(animals, 'stable')'
  pilotIDJPNormProfile('Animal', animal, 'NBootstrap', opts.NBootstrap, 'Verbose', opts.Verbose);
end
end
