function T = plotIncDecScaleCorrelation(kernelDir, varargin)
% plotIncDecScaleCorrelation
%
% Extract DEC and INC scale values across kernel files and test whether
% they are correlated across sessions.
%
% Assumes:
%   compStats.scale is [nSideTypes x 2], with columns:
%       1 = DEC
%       2 = INC
%   sideTypeNames gives the row labels
%   header.probeDirDeg contains the probe direction offset in degrees
%
% Example:
%   T = plotIncDecScaleCorrelation( ...
%       '/Users/Shared/Data/IDReadout/Data/Kernels', ...
%       'ScaleSideType', 'change', ...
%       'ProbeOffsetDeg', 45);

p = inputParser;
addParameter(p, 'ScaleSideType', 'change', @(x)ischar(x)||isstring(x));
addParameter(p, 'ProbeOffsetDeg', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'UseAbsOffset', true, @islogical);
addParameter(p, 'MakePlot', true, @islogical);
parse(p, varargin{:});

scaleSideType = char(p.Results.ScaleSideType);
probeOffsetDeg = p.Results.ProbeOffsetDeg;
useAbsOffset   = p.Results.UseAbsOffset;
makePlot       = p.Results.MakePlot;

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

    probeDir = localExtractScalar(S.header.probeDirDeg);
    if isnan(probeDir)
        fprintf('Skipping %s: could not read probeDirDeg\n', files(iFile).name);
        continue;
    end

    if isempty(probeOffsetDeg)
        keepFile = true;   % pool across all probe directions
    elseif useAbsOffset
        keepFile = abs(probeDir) == probeOffsetDeg;
    else
        keepFile = probeDir == probeOffsetDeg;
    end
    
    if ~keepFile
        continue;
    end

    sideIdx = find(strcmpi(S.sideTypeNames, scaleSideType), 1);
    if isempty(sideIdx)
        fprintf('Skipping %s: side type "%s" not found\n', files(iFile).name, scaleSideType);
        continue;
    end

    scaleMat = S.compStats.scale;
    if size(scaleMat,2) < 2
        fprintf('Skipping %s: compStats.scale has wrong shape\n', files(iFile).name);
        continue;
    end

    decScale = scaleMat(sideIdx, 1);
    incScale = scaleMat(sideIdx, 2);

    datFile = "";
    if isfield(S.header, 'fileName')
        datFile = string(localExtractString(S.header.fileName));
    end

    rows(end+1,:) = {files(iFile).name, datFile, probeDir, decScale, incScale}; %#ok<AGROW>
end

if isempty(rows)
    warning('No matching sessions found.');
    T = table;
    return;
end

T = cell2table(rows, 'VariableNames', ...
    {'kernelFile','datFile','probeDirDeg','decScale','incScale'});

keep = isfinite(T.decScale) & isfinite(T.incScale);
T = T(keep,:);

disp(T);

n = height(T);
if isempty(probeOffsetDeg)
    fprintf('\nFound %d sessions pooled across all probe directions\n', n);
else
    fprintf('\nFound %d sessions for probeDirDeg = %g\n', n, probeOffsetDeg);
end

if n >= 3
    [rP, pP] = corr(T.decScale, T.incScale, 'Type', 'Pearson');
    [rS, pS] = corr(T.decScale, T.incScale, 'Type', 'Spearman');

    fprintf('Pearson:  r = %.4f, p = %.4g\n', rP, pP);
    fprintf('Spearman: r = %.4f, p = %.4g\n', rS, pS);
else
    rP = NaN; pP = NaN; rS = NaN; pS = NaN;
    fprintf('Too few sessions for a meaningful correlation test.\n');
end

if makePlot
    figure(301); clf; hold on;
    plot(T.decScale, T.incScale, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7);

    for i = 1:n
        text(T.decScale(i), T.incScale(i), sprintf('  %d', i), ...
            'FontSize', 9, 'VerticalAlignment', 'middle');
    end

    lo = min([T.decScale; T.incScale]);
    hi = max([T.decScale; T.incScale]);
    pad = max(0.03, 0.05*(hi-lo));
    xlim([lo-pad hi+pad]);
    ylim([lo-pad hi+pad]);
    plot([lo-pad hi+pad], [lo-pad hi+pad], 'k--');

    xlabel('DEC scale');
    ylabel('INC scale');
    if isempty(probeOffsetDeg)
      title(sprintf('%s scales pooled across all probe directions', scaleSideType));
    else
      title(sprintf('%s scales across %g^\\circ probe sessions', scaleSideType, probeOffsetDeg));
    end
    axis square;
    box off;

    if n >= 3
        txt = sprintf('n = %d\nPearson r = %.3f (p = %.3g)\nSpearman r = %.3f (p = %.3g)', ...
            n, rP, pP, rS, pS);
    else
        txt = sprintf('n = %d', n);
    end

    xl = xlim; yl = ylim;
    text(xl(1)+0.05*range(xl), yl(2)-0.05*range(yl), txt, ...
        'VerticalAlignment', 'top', 'FontSize', 10);
end
end

function x = localExtractScalar(v)
% Handle either plain numeric or Lablib-style structs with .data
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