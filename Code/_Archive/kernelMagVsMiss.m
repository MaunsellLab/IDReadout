function kernelMagVsMiss(baseFolder)

% excludedFiles = {'IDReadout_Meetz_20260113', 'IDReadout_Meetz_20260114', ...
%   'IDReadout_Meetz_20260114_2', 'IDReadout_Meetz_20260114_3', ...
%   'IDReadout_Meetz_20260114_4', 'IDReadout_Meetz_20260115', ...
%   'IDReadout_Meetz_20260116', 'IDReadout_Meetz_20260209', ...
%   'IDReadout_Meetz_20260210'};
excludedFiles = {};

if nargin < 1 || isempty(baseFolder)
    baseFolder = folderPath();
end

dataFolder = baseFolder + "/Kernels/";
if ~exist(dataFolder, 'dir')
    error('Kernels folder not found: %s', dataFolder);
end

matFiles = dir(fullfile(dataFolder, '*.mat'));

% ---------- labels ----------
% Adjust these if your sideType / pref-opp labels differ.
sideLabels = {'Change', 'NoChange', 'Diff', 'RF', 'Opp'};
stepLabels = {'INC', 'DEC'};
polLabels  = {'Pref', 'Anti'};   % rename if needed

% ---------- accumulators ----------
fileNames = {};
sessionDates = {};
kIntAll = [];        % [nSessions x 5 x 2 x 2]
nMissAll = [];       % [nSessions x 2]
nRFMissAll = [];     % [nSessions x 2]
nOppMissAll = [];    % [nSessions x 2]

nSessions = 0;

for f = 1:length(matFiles)
    fileName = matFiles(f).name;
    [~, baseName, ~] = fileparts(fileName);

    if endsWith(fileName, '_fileInfo.mat')
        continue;
    end
    if ~isempty(excludedFiles) && ismember(baseName, excludedFiles)
        continue;
    end

    S = load(dataFolder + fileName);

    if ~isfield(S, 'header') || ~isfield(S, 'compStats') || ~isfield(S, 'hitStats')
        fprintf('Skipping   %s -- missing required variables\n', fileName);
        continue;
    end

    header = S.header;

    % Keep same filtering used in kernelAverage
    if header.prefNoiseCohPC.data ~= 10
        fprintf("Skipping   %s (%d of %d) -- prefNoiseCohPC is %.0f\n", ...
            fileName, f, length(matFiles), header.prefNoiseCohPC.data);
        continue;
    end

    fprintf("Processing %s (%d of %d)\n", fileName, f, length(matFiles));

    compStats = S.compStats;
    hitStats  = S.hitStats;

    if ~isfield(compStats, 'kIntegrals')
        fprintf('Skipping   %s -- compStats.kIntegrals missing\n', fileName);
        continue;
    end
    if ~isfield(hitStats, 'nTrials') || ~isfield(hitStats, 'nHits') || ...
       ~isfield(hitStats, 'nRFTrials') || ~isfield(hitStats, 'nRFHits')
        fprintf('Skipping   %s -- hitStats missing required fields\n', fileName);
        continue;
    end

    nSessions = nSessions + 1;

    kIntAll(nSessions, :, :, :) = compStats.kIntegrals;

    % Overall misses by INC / DEC block
    nMissAll(nSessions, :) = hitStats.nTrials - hitStats.nHits;

    % RF misses by INC / DEC block
    nRFMissAll(nSessions, :) = hitStats.nRFTrials - hitStats.nRFHits;

    % Opp misses by INC / DEC block
    oppTrials = hitStats.nTrials - hitStats.nRFTrials;
    oppHits   = hitStats.nHits   - hitStats.nRFHits;
    nOppMissAll(nSessions, :) = oppTrials - oppHits;

    fileNames{nSessions,1} = fileName;

    if isfield(header, 'date') && isfield(header.date, 'data')
        sessionDates{nSessions,1} = header.date.data;
    else
        sessionDates{nSessions,1} = '';
    end
end

if nSessions == 0
    error('No usable session kernel files found.');
end

% -------------------------------------------------------------------------
% Choose which miss count to use on x-axis:
%   nMissAll   = total misses in each INC/DEC block
%   nRFMissAll = misses on RF-side trials
%   nOppMissAll = misses on Opp-side trials
%
% For your current diagnostic, total misses is probably the most natural.
% -------------------------------------------------------------------------
xMiss = nMissAll;

% Use absolute integral as kernel magnitude.
yMag = abs(kIntAll);

% If you want signed integrals instead, replace the line above with:
% yMag = kIntAll;

figure(4); clf;
set(gcf, 'Name', 'Kernel Magnitude vs Miss Count');

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for s = 1:2          % INC / DEC
    for p = 1:2      % Pref / Anti
        nexttile;
        hold on;

        for sideType = 1:5
            x = xMiss(:, s);
            y = squeeze(yMag(:, sideType, s, p));

            good = isfinite(x) & isfinite(y);
            if any(good)
                scatter(x(good), y(good), 45, 'filled', ...
                    'DisplayName', sideLabels{sideType});

                % Optional least-squares fit line
                if sum(good) >= 2
                    pp = polyfit(x(good), y(good), 1);
                    xx = linspace(min(x(good)), max(x(good)), 100);
                    yy = polyval(pp, xx);
                    plot(xx, yy, 'LineWidth', 1.0, ...
                        'HandleVisibility', 'off');
                end
            end
        end

        xlabel('Miss trials');
        ylabel('|Kernel integral|');
        title(sprintf('%s  %s', stepLabels{s}, polLabels{p}));
        box off;
        grid on;
    end
end

lgd = legend('Location', 'bestoutside');
lgd.Box = 'off';

sgtitle(sprintf('Kernel Magnitude vs Miss Count (%d sessions)', nSessions));

end