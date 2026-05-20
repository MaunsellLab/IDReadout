function summary = checkPrefProbeNoiseCorrelation(dataRoot)
% checkPrefProbeNoiseCorrelation
%
% Standalone sanity check for independence of preferred and probe coherence
% noise streams in the MT readout experiment.
%
% Walks:
%     Data/probe*/NoiseMatrices/*.mat
%
% Each file is expected to contain:
%     prefNoiseByPatch   [2 x nFrames x nTrials]
%     probeNoiseByPatch  [2 x nFrames x nTrials]
%
% For each file, the function vectorizes both matrices and computes one
% Pearson correlation coefficient. It then reports the mean and SEM across
% files, along with a rough expected SEM assuming independent binary noise.
%
% Usage:
%     summary = checkPrefProbeNoiseCorrelation();
%     summary = checkPrefProbeNoiseCorrelation('/path/to/Data');

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(folderPath(), 'Data');
end

assert(isfolder(dataRoot), 'Data root does not exist: %s', dataRoot);

probeDirs = dir(fullfile(dataRoot, 'probe*'));
probeDirs = probeDirs([probeDirs.isdir]);

rows = struct( ...
    'fileName', {}, ...
    'probeTag', {}, ...
    'nObs', {}, ...
    'r', {}, ...
    'prefMean', {}, ...
    'probeMean', {}, ...
    'prefSD', {}, ...
    'probeSD', {} );

for p = 1:numel(probeDirs)
    probeTag = probeDirs(p).name;
    noiseDir = fullfile(dataRoot, probeTag, 'NoiseMatrices');

    if ~isfolder(noiseDir)
        continue
    end

    files = dir(fullfile(noiseDir, '*.mat'));

    for f = 1:numel(files)
        filePath = fullfile(noiseDir, files(f).name);

        S = load(filePath, 'prefNoiseByPatch', 'probeNoiseByPatch');

        if ~isfield(S, 'prefNoiseByPatch') || ~isfield(S, 'probeNoiseByPatch')
            warning('Skipping %s: missing prefNoiseByPatch or probeNoiseByPatch.', filePath);
            continue
        end

        pref  = S.prefNoiseByPatch;
        probe = S.probeNoiseByPatch;

        assert(isequal(size(pref), size(probe)), ...
            'Size mismatch in %s: pref is %s, probe is %s.', ...
            filePath, mat2str(size(pref)), mat2str(size(probe)));

        assert(ndims(pref) == 3 && size(pref,1) == 2, ...
            'Unexpected prefNoiseByPatch size in %s: %s.', ...
            filePath, mat2str(size(pref)));

        x = pref(:);
        y = probe(:);

        valid = isfinite(x) & isfinite(y);
        x = x(valid);
        y = y(valid);

        if numel(x) < 4
            warning('Skipping %s: fewer than 4 valid observations.', filePath);
            continue
        end

        if std(x) == 0 || std(y) == 0
            warning('Skipping %s: zero variance in pref or probe noise.', filePath);
            continue
        end

        C = corrcoef(x, y);
        r = C(1,2);

        k = numel(rows) + 1;
        rows(k).fileName  = filePath;
        rows(k).probeTag  = probeTag;
        rows(k).nObs      = numel(x);
        rows(k).r         = r;
        rows(k).prefMean  = mean(x);
        rows(k).probeMean = mean(y);
        rows(k).prefSD    = std(x);
        rows(k).probeSD   = std(y);
    end
end

assert(~isempty(rows), 'No usable NoiseMatrices files found under %s.', dataRoot);

T = struct2table(rows);

rVals = T.r;
nFiles = height(T);

meanR = mean(rVals);
semR  = std(rVals, 0) / sqrt(nFiles);

% Rough null expectation:
% For independent streams, the sampling SD of a Pearson r is approximately
% 1/sqrt(N - 3) using the Fisher-z approximation. This is essentially the
% same as 1/sqrt(N) for these large N values.
expectedVarPerFile = 1 ./ max(T.nObs - 3, 1);

% Expected SEM for the unweighted mean of the per-file r values.
expectedSEMMean = sqrt(sum(expectedVarPerFile)) / nFiles;

% Also report an intuitive typical single-file null SD.
typicalExpectedFileSD = mean(sqrt(expectedVarPerFile));

fprintf('\nPreferred/probe noise correlation sanity check\n');
fprintf('Data root: %s\n', dataRoot);
fprintf('Files analyzed: %d\n', nFiles);
fprintf('Mean r across files: %.6g\n', meanR);
fprintf('SEM across files:    %.6g\n', semR);
fprintf('Expected SEM under independence, approx: %.6g\n', expectedSEMMean);
fprintf('Typical expected SD of one file''s r, approx: %.6g\n', typicalExpectedFileSD);
fprintf('Mean observations per file: %.1f\n', mean(T.nObs));
fprintf('Range observations per file: %d to %d\n', min(T.nObs), max(T.nObs));

fprintf('\nLargest absolute correlations:\n');
[~, ix] = sort(abs(T.r), 'descend');
disp(T(ix(1:min(10,nFiles)), {'probeTag','nObs','r','prefMean','probeMean','fileName'}));

summary = struct();
summary.dataRoot = dataRoot;
summary.fileTable = T;
summary.nFiles = nFiles;
summary.meanR = meanR;
summary.semR = semR;
summary.expectedSEMMean = expectedSEMMean;
summary.typicalExpectedFileSD = typicalExpectedFileSD;

end