function T = plotScaleVsPrefDir(kernelDir, varargin)
% plotScaleVsPrefDir
%
% Plot per-session scale values against session probe/drift direction
% (taken from header.prefDirDeg).
%
% Defaults:
%   - step type: INC
%   - side type: change
%
% Example:
%   T = plotScaleVsPrefDir( ...
%       '/Users/Shared/Data/IDReadout/Data/Kernels');
%
%   T = plotScaleVsPrefDir( ...
%       '/Users/Shared/Data/IDReadout/Data/Kernels', ...
%       'ScaleSideType', 'diff', ...
%       'StepType', 'dec');

p = inputParser;
addParameter(p, 'ScaleSideType', 'change', @(x)ischar(x)||isstring(x));
addParameter(p, 'StepType', 'inc', @(x)ischar(x)||isstring(x));
addParameter(p, 'MakePlot', true, @islogical);
addParameter(p, 'JitterDeg', 2.5, @isscalar);
parse(p, varargin{:});

scaleSideType = char(p.Results.ScaleSideType);
stepType      = lower(char(p.Results.StepType));
makePlot      = p.Results.MakePlot;
jitterDeg     = p.Results.JitterDeg;

switch stepType
    case 'dec'
        stepIdx = 1;
    case 'inc'
        stepIdx = 2;
    otherwise
        error('StepType must be ''dec'' or ''inc''.');
end

files = dir(fullfile(kernelDir, '*.mat'));
if isempty(files)
    error('No kernel files found in %s', kernelDir);
end

rows = {};

for iFile = 1:numel(files)
    f = fullfile(files(iFile).folder, files(iFile).name);
    S = load(f, 'header', 'compStats', 'sideTypeNames');

    if ~isfield(S, 'header') || ~isfield(S, 'compStats') || ~isfield(S, 'sideTypeNames')
        fprintf('Skipping %s: missing needed variables\n', files(iFile).name);
        continue;
    end

    if ~isfield(S.compStats, 'scale')
        fprintf('Skipping %s: compStats.scale missing\n', files(iFile).name);
        continue;
    end

    prefDir = localExtractScalar(S.header.prefDirDeg);
    if isnan(prefDir)
        fprintf('Skipping %s: could not read prefDirDeg\n', files(iFile).name);
        continue;
    end

    sideIdx = find(strcmpi(S.sideTypeNames, scaleSideType), 1);
    if isempty(sideIdx)
        fprintf('Skipping %s: side type "%s" not found\n', files(iFile).name, scaleSideType);
        continue;
    end

    scaleMat = S.compStats.scale;
    if size(scaleMat,2) < stepIdx
        fprintf('Skipping %s: compStats.scale has wrong shape\n', files(iFile).name);
        continue;
    end

    scaleVal = scaleMat(sideIdx, stepIdx);

    datFile = "";
    if isfield(S.header, 'fileName')
        datFile = string(localExtractString(S.header.fileName));
    end

    rows(end+1,:) = {files(iFile).name, datFile, prefDir, scaleVal}; %#ok<AGROW>
end

if isempty(rows)
    warning('No usable sessions found.');
    T = table;
    return;
end

T = cell2table(rows, 'VariableNames', ...
    {'kernelFile','datFile','prefDirDeg','scale'});

T = T(isfinite(T.prefDirDeg) & isfinite(T.scale), :);

fprintf('\nFound %d sessions\n', height(T));
disp(T(:, {'prefDirDeg','scale','datFile'}));

if makePlot
    figure(302); clf; hold on;

    % jittered individual sessions
    x = T.prefDirDeg + (rand(height(T),1)-0.5)*2*jitterDeg;
    plot(x, T.scale, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    % mean by 45-degree bins
    dirWrapped = mod(T.prefDirDeg, 360);
    binCtrs = 0:45:315;

    meanScale = nan(size(binCtrs));
    semScale  = nan(size(binCtrs));
    nPerBin   = zeros(size(binCtrs));

    % Assign each session to nearest 45-degree bin center
    % Equivalent to bins spanning center +/- 22.5 deg
    binIdx = mod(round(dirWrapped / 45), 8) + 1;

    for i = 1:numel(binCtrs)
        vals = T.scale(binIdx == i);
        nPerBin(i) = numel(vals);

        if ~isempty(vals)
            meanScale(i) = mean(vals);
            if numel(vals) > 1
                semScale(i) = std(vals) / sqrt(numel(vals));
            end
        end
    end

    errorbar(binCtrs, meanScale, semScale, 'ks-', ...
        'LineWidth', 1.5, 'MarkerFaceColor', 'w', 'MarkerSize', 7);    xlabel('Drift direction (deg)');

    xlim([-10 325]);
    xticks(0:45:315);
    ylabel(sprintf('%s scale', upper(stepType)));
    title(sprintf('%s, %s scale vs direct direction', scaleSideType, upper(stepType)));
    box off;

    % optional horizontal zero line
    xl = xlim;
    plot(xl, [0 0], 'k:');
    xlim(xl);
end
end

function x = localExtractScalar(v)
x = NaN;
if isnumeric(v) && isscalar(v)
    x = double(v);
elseif isstruct(v) && isfield(v, 'data')
    if isnumeric(v.data) && isscalar(v.data)
        x = double(v.data);
    end
end
end

function s = localExtractString(v)
s = "";
if ischar(v) || isstring(v)
    s = string(v);
elseif isstruct(v) && isfield(v, 'data')
    s = string(v.data);
end
end