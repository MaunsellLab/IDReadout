function S = summarizeStrategyModels(basePath)

% ---- basePath handling ----
if nargin < 1 || isempty(basePath)
    basePath = folderPath();
end
if isstring(basePath)
    basePath = char(basePath);
end

modelFolder   = fullfile(basePath,'Models');
summaryFolder = fullfile(basePath,'Summaries');

if ~exist(summaryFolder,'dir')
    mkdir(summaryFolder);
end

% ---- discover model files ----
d = dir(fullfile(modelFolder,'*.mat'));

if isempty(d)
    fprintf('No model files found in %s\n',modelFolder);
    S = [];
    return
end

modelFiles = fullfile({d.folder},{d.name});

% ============================================================
% Build session summary table
% ============================================================

rows = struct([]);

for i = 1:numel(modelFiles)

    try
        tmp = load(modelFiles{i},'M');
        row = tmp.M.summaryRow;

        rows = [rows; row]; %#ok<AGROW>

    catch ME
        warning('Could not read %s (%s)',modelFiles{i},ME.message)
    end

end
if isempty(rows)
    warning('summarizeStrategyModels:NoValidRows', ...
        'No valid summary rows were loaded from %s', modelFolder);
    S = [];
    return
end
sessionTable = struct2table(rows);

% ---- sort by session name for stability ----
if ismember('session',sessionTable.Properties.VariableNames)
    sessionTable = sortrows(sessionTable,'session');
end

% ============================================================
% pooled analysis
% ============================================================

[tblPool, info] = collectStrategyTables(modelFiles);

P = poolStrategyForensics(tblPool,'verbose',false);

% ============================================================
% Save outputs
% ============================================================

csvFile = fullfile(summaryFolder,'strategy_session_summary.csv');
matFile = fullfile(summaryFolder,'strategy_session_summary.mat');
poolFile = fullfile(summaryFolder,'strategy_pooled_summary.mat');

writetable(sessionTable,csvFile)

save(matFile,'sessionTable','info')
save(poolFile,'P')

% ============================================================
% short console summary
% ============================================================

fprintf('\nStrategy summary\n')
fprintf('Sessions: %d\n',height(sessionTable))
fprintf('Trials pooled: %d\n',height(tblPool))

if ismember('pCorrect',sessionTable.Properties.VariableNames)
    fprintf('Mean session P(correct): %.3f\n',mean(sessionTable.pCorrect,'omitnan'))
end

if isfield(P,'summary')
    fprintf('\nPooled model deviance explained:\n')
    fprintf('  m1: %.5f\n',P.summary.devExplained_p1)
    fprintf('  m2: %.5f\n',P.summary.devExplained_p2)
end

% ============================================================
% return struct
% ============================================================

S = struct();
S.sessionTable = sessionTable;
S.pooled = P;
S.info = info;

plotStrategySummary(S, basePath);

end