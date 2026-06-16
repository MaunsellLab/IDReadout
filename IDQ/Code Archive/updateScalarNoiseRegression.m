function batch = updateScalarNoiseRegression(replace, varargin)
% updateScalarNoiseRegression  Build scalar noise regression analyses.
%
%   batch = updateScalarNoiseRegression(replace)
%
% Runs low-dimensional logistic regressions on probe-session files:
%
%   Data/ProbeXX/ProbeSessions/*.mat
%
% and writes outputs to:
%
%   Data/ProbeXX/Regression/*_scalarNoiseRegression.mat
%   Plots/Probes/ProbeXX/Regression/*_scalarNoiseRegression.pdf
%
% The primary analysis uses mean preferred-direction coherence noise over the
% 0-250 ms post-step window, represented by the last nStepFrames entries in
% prefNoiseByPatch. Because prefNoiseByPatch stores the realized noise value
% on each video frame, noise-frame phase is already handled by averaging over
% video frames.
%
% Inputs
%   replace : logical. If false, skip files with existing .mat and .pdf.
%
% Name-value parameters
%   'WinStartMS' : default 0. Currently only supports 0 ms.
%   'WinDurMS'   : default 250.
%   'NBins'      : default 10.
%
% Output
%   batch : summary struct saved to Data/scalarNoiseRegression_batchSummary.mat

% cleanupObj = initProjectPath(); %#ok<NASGU>
if nargin < 1 || isempty(replace)
    replace = false;
end

p = inputParser;
p.addParameter('WinStartMS', 0, @isscalar);
p.addParameter('WinDurMS', 250, @isscalar);
p.addParameter('NBins', 10, @isscalar);
p.parse(varargin{:});

pathRoot = domainFolder(mfilename('fullpath'));
dataRoot = fullfile(pathRoot, 'Data');
plotRoot = fullfile(pathRoot, 'Plots', 'Probes');

probeDirs = dir(fullfile(dataRoot, 'Probe*'));
probeDirs = probeDirs([probeDirs.isdir]);

batch = struct();
batch.params = p.Results;
batch.replace = replace;
batch.files = {};
batch.regPaths = {};
batch.plotPaths = {};
batch.ok = [];
batch.messages = {};

if isempty(probeDirs)
    fprintf('No Probe* folders found in %s\n', dataRoot);
    return;
end

for iProbe = 1:numel(probeDirs)

    probeTag = probeDirs(iProbe).name;
    probeDataFolder = fullfile(dataRoot, probeTag);
    probeSessionFolder = fullfile(probeDataFolder, 'ProbeSessions');

    if ~isfolder(probeSessionFolder)
        continue;
    end

    regressionFolder = validFolder(fullfile(probeDataFolder, 'Regression'));
    plotFolder = validFolder(fullfile(plotRoot, probeTag, 'Regression'));

    sessionFiles = dir(fullfile(probeSessionFolder, '*.mat'));

    for iFile = 1:numel(sessionFiles)

        probeSessionPath = fullfile(sessionFiles(iFile).folder, sessionFiles(iFile).name);
        [~, baseName] = fileparts(probeSessionPath);

        regPath = fullfile(regressionFolder, sprintf('%s_scalarNoiseRegression.mat', baseName));
        plotPath = fullfile(plotFolder, sprintf('%s_scalarNoiseRegression.pdf', baseName));

        batch.files{end+1, 1} = probeSessionPath; %#ok<AGROW>
        batch.regPaths{end+1, 1} = regPath; %#ok<AGROW>
        batch.plotPaths{end+1, 1} = plotPath; %#ok<AGROW>
        batch.ok(end+1, 1) = false; %#ok<AGROW>
        batch.messages{end+1, 1} = ''; %#ok<AGROW>
        iBatch = numel(batch.files);

        if ~replace && isfile(regPath) && isfile(plotPath)
            fprintf('Skipping existing regression: %s [%s]\n', baseName, probeTag);
            batch.ok(iBatch) = true;
            batch.messages{iBatch} = 'existing outputs; skipped';
            continue;
        end

        fprintf('Processing scalar noise regression: %s [%s]\n', baseName, probeTag);

        try
            reg = computeScalarNoiseRegression(probeSessionPath, ...
                'WinStartMS', p.Results.WinStartMS, ...
                'WinDurMS', p.Results.WinDurMS, ...
                'NBins', p.Results.NBins);

            save(regPath, 'reg', '-v7.3');

            fig = plotScalarNoiseRegression(reg);
            exportgraphics(fig, plotPath, 'ContentType', 'vector');
            close(fig);

            batch.ok(iBatch) = true;
            batch.messages{iBatch} = 'ok';

            fprintf('  beta change pref:   %.4g\n', reg.effectSummary.bothPref.DchangePref.beta);
            fprintf('  beta no-change pref: %.4g\n', reg.effectSummary.bothPref.DnoChangePref.beta);

        catch ME
            warning('Failed on %s:\n%s', probeSessionPath, ME.message);
            batch.messages{iBatch} = ME.message;
        end
    end
end

batchSummaryPath = fullfile(dataRoot, 'scalarNoiseRegression_batchSummary.mat');
save(batchSummaryPath, 'batch');

fprintf('\nScalar noise regression complete: %d/%d successful.\n', sum(batch.ok), numel(batch.ok));
fprintf('Batch summary saved: %s\n', batchSummaryPath);

end
