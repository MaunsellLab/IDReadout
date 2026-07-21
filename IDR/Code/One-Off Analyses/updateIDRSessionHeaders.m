function updateIDRSessionHeaders(varargin)
% updateIDRSessionHeaders  Add compact behavior tables to converted IDR files.
%
% updateIDRSessionHeaders(...)
%
% Selects parent FullSessions with selectAnalysisFiles, loads header, trials,
% and the existing sessionHeader, rebuilds sessionHeader with the standalone
% makeSessionHeader function, and appends only the revised sessionHeader to
% each MAT-file.
%
% Name-value arguments:
%   Animal          Animal name or 'All' (default)
%   Replace         Rebuild even when behaviorTrialTable already exists
%                   (default false)
%   ReportExcluded  Forwarded to selectAnalysisFiles (default false)
%
% Example:
%   updateIDRSessionHeaders('Animal', 'Neesha');
%
% The full trials cell array is loaded once per file during this update.
% Subsequent psychometric inventories can use sessionHeader.behaviorTrialTable
% without loading trials.

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Animal', 'All', @(x) ischar(x) || isstring(x));
addParameter(p, 'Replace', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'ReportExcluded', false, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opts = p.Results;

domainPath = domainFolder(mfilename('fullpath'));
fullSessionFolder = fullfile(domainPath, 'Data', 'FullSessions');

[filePaths, fileInfo] = selectAnalysisFiles( ...
  fullSessionFolder, ...
  'Animal', opts.Animal, ...
  'ReportExcluded', opts.ReportExcluded);

if isempty(filePaths)
  fprintf('updateIDRSessionHeaders: no FullSession files selected.\n');
  return;
end

nUpdated = 0;
nSkipped = 0;

for iFile = 1:numel(filePaths)
  filePath = filePaths{iFile};

  % Check the lightweight header before loading the large trials cell array.
  S = load(filePath, 'sessionHeader');
  if isfield(S, 'sessionHeader') && ...
      isfield(S.sessionHeader, 'behaviorTrialTable') && ...
      ~opts.Replace
    fprintf('Skipping %s: behaviorTrialTable already exists.\n', ...
      fileInfo.fileName{iFile});
    nSkipped = nSkipped + 1;
    continue;
  end

  fprintf('Updating %s ...\n', fileInfo.fileName{iFile});

  S = load(filePath, 'header', 'trials', 'sessionHeader');

  if ~isfield(S, 'header') || ~isfield(S, 'trials')
    error('updateIDRSessionHeaders:MissingVariables', ...
      'File %s must contain header and trials.', filePath);
  end

  % The existing sessionHeader supplies the session-level trialMeta fields
  % produced by convertIDRData. If it is absent, makeSessionHeader will use
  % defaults, but converted FullSessions are expected to contain it.
  if isfield(S, 'sessionHeader')
    trialMeta = S.sessionHeader;
  else
    trialMeta = struct();
  end

  sessionHeader = makeSessionHeader(S.header, trialMeta, S.trials); %#ok<NASGU>

  % Append only the revised header. This avoids rewriting the large file.
  save(filePath, 'sessionHeader', '-append');

  nUpdated = nUpdated + 1;
end

fprintf(['updateIDRSessionHeaders complete: %d updated, %d skipped, ' ...
  '%d selected.\n'], nUpdated, nSkipped, numel(filePaths));
end
