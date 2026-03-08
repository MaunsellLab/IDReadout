function [tblAll, info] = collectStrategyTables(modelFiles, varargin)
% collectStrategyTables  Concatenate per-session strategy tables from Models/.
%
% Usage:
%   tblAll = collectStrategyTables(modelFiles)
%   [tblAll, info] = collectStrategyTables(modelFiles)
%
% Input:
%   modelFiles   cell array of full paths to per-session model .mat files
%                Each file must contain struct M with field M.tblAll.
%
% Optional name/value:
%   'verbose'    logical, default true
%
% Output:
%   tblAll       concatenated trial table across sessions
%   info         struct with bookkeeping information
%
% Notes:
%   - Adds the following variables to each session table if absent:
%         sessionID      categorical session label
%         sessionIndex   numeric session index in input list
%         sessionFile    string session filename
%         modelFile      string full path to model file
%   - Rows from files that fail validation are skipped with warning.

% ----------------------------
% parse inputs
% ----------------------------
p = inputParser;
p.addRequired('modelFiles', @(x) iscell(x) || isstring(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(modelFiles, varargin{:});

verbose = p.Results.verbose;

if isstring(modelFiles)
    modelFiles = cellstr(modelFiles);
end

nFiles = numel(modelFiles);

% ----------------------------
% initialize outputs
% ----------------------------
tblAll = table();

info = struct();
info.modelFilesRequested = modelFiles(:);
info.loaded = false(nFiles,1);
info.skipped = false(nFiles,1);
info.skipReason = strings(nFiles,1);
info.nTrialsPerSession = nan(nFiles,1);
info.sessionLabels = strings(nFiles,1);

if nFiles == 0
    if verbose
        fprintf('collectStrategyTables: no model files supplied.\n');
    end
    return;
end

% ----------------------------
% loop over files
% ----------------------------
for iFile = 1:nFiles
    thisFile = modelFiles{iFile};

    try
        S = load(thisFile, 'M');
    catch ME
        warning('collectStrategyTables:LoadFailed', ...
            'Could not load "%s": %s', thisFile, ME.message);
        info.skipped(iFile) = true;
        info.skipReason(iFile) = "load failed";
        continue;
    end

    if ~isfield(S, 'M') || ~isstruct(S.M)
        warning('collectStrategyTables:MissingM', ...
            'File "%s" does not contain struct M.', thisFile);
        info.skipped(iFile) = true;
        info.skipReason(iFile) = "missing struct M";
        continue;
    end

    M = S.M;

    if ~isfield(M, 'tblAll') || ~istable(M.tblAll) || isempty(M.tblAll)
        warning('collectStrategyTables:MissingTblAll', ...
            'File "%s" has no valid M.tblAll.', thisFile);
        info.skipped(iFile) = true;
        info.skipReason(iFile) = "missing or empty M.tblAll";
        continue;
    end

    T = M.tblAll;

    % session label preference:
    %   1) M.sessionFile if available
    %   2) basename of model file
    if isfield(M, 'sessionFile') && ~isempty(M.sessionFile)
        sessionLabel = string(M.sessionFile);
    else
        [~, base, ~] = fileparts(thisFile);
        sessionLabel = string(base);
    end

    info.sessionLabels(iFile) = sessionLabel;
    info.nTrialsPerSession(iFile) = height(T);

    % append bookkeeping variables
    nRows = height(T);

    if ~ismember('sessionID', T.Properties.VariableNames)
        T.sessionID = categorical(repmat(sessionLabel, nRows, 1));
    else
        T.sessionID = categorical(string(T.sessionID));
    end

    if ~ismember('sessionIndex', T.Properties.VariableNames)
        T.sessionIndex = repmat(iFile, nRows, 1);
    end

    if ~ismember('sessionFile', T.Properties.VariableNames)
        if isfield(M, 'sessionFile') && ~isempty(M.sessionFile)
            T.sessionFile = repmat(string(M.sessionFile), nRows, 1);
        else
            T.sessionFile = repmat(sessionLabel, nRows, 1);
        end
    else
        T.sessionFile = string(T.sessionFile);
    end

    if ~ismember('modelFile', T.Properties.VariableNames)
        T.modelFile = repmat(string(thisFile), nRows, 1);
    else
        T.modelFile = string(T.modelFile);
    end

    % concatenate
    if isempty(tblAll)
        tblAll = T;
    else
        try
            tblAll = [tblAll; T]; %#ok<AGROW>
        catch ME
            warning('collectStrategyTables:ConcatFailed', ...
                ['Could not concatenate table from "%s".\n' ...
                 'This usually means variable names/types differ across sessions.\n' ...
                 'Error: %s'], thisFile, ME.message);
            info.skipped(iFile) = true;
            info.skipReason(iFile) = "concatenation failed";
            continue;
        end
    end

    info.loaded(iFile) = true;
end

% ----------------------------
% final info
% ----------------------------
info.nFilesRequested = nFiles;
info.nFilesLoaded = sum(info.loaded);
info.nFilesSkipped = sum(info.skipped);

if isempty(tblAll)
    info.nTrialsTotal = 0;
    info.nSessionsLoaded = 0;
else
    info.nTrialsTotal = height(tblAll);
    info.nSessionsLoaded = numel(categories(tblAll.sessionID));
end

if verbose
    fprintf('collectStrategyTables:\n');
    fprintf('  requested files: %d\n', info.nFilesRequested);
    fprintf('  loaded files:    %d\n', info.nFilesLoaded);
    fprintf('  skipped files:   %d\n', info.nFilesSkipped);
    fprintf('  pooled sessions: %d\n', info.nSessionsLoaded);
    fprintf('  pooled trials:   %d\n', info.nTrialsTotal);
end

end