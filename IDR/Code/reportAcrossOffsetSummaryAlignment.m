function reportAcrossOffsetSummaryAlignment()
% reportAcrossOffsetSummaryAlignment
%
% Reports whether the across-offset beta and kernel summaries contain the
% same probe-session files.
%
% Expected files:
%   IDR/Data/AcrossOffsetSummaries/IDR_BetaSummary.mat
%   IDR/Data/AcrossOffsetSummaries/IDR_KernelSummary.mat
%
% Beta summary file:
%   betaSummary.sessionRecords(:).fileName
%
% Kernel summary file:
%   acrossOffsetSummary.fileInfo.fileName
%
% Beta file names are normalized by removing:
%   _scalarNoiseRegression
%
% Comparisons are case-insensitive.

    projectRoot = currentProject().RootFolder;
    summaryFolder = fullfile(projectRoot, 'IDR', 'Data', 'AcrossOffsetSummaries');
    betaSummaryPath = fullfile(summaryFolder, 'IDR_BetaSummary.mat');
    kernelSummaryPath = fullfile(summaryFolder, 'IDR_KernelSummary.mat');

    fprintf('\nAcross-offset summary alignment report\n');
    fprintf('--------------------------------------\n');

    if ~isfile(betaSummaryPath)
        error('Beta summary file not found:\n  %s', betaSummaryPath);
    end

    if ~isfile(kernelSummaryPath)
        error('Kernel summary file not found:\n  %s', kernelSummaryPath);
    end

    betaData = load(betaSummaryPath, 'betaSummary');
    kernelData = load(kernelSummaryPath, 'acrossOffsetSummary');

    if ~isfield(betaData, 'betaSummary')
        error('File does not contain betaSummary:\n  %s', betaSummaryPath);
    end

    if ~isfield(kernelData, 'acrossOffsetSummary')
        error('File does not contain acrossOffsetSummary:\n  %s', kernelSummaryPath);
    end

    betaSummary = betaData.betaSummary;
    acrossOffsetSummary = kernelData.acrossOffsetSummary;

    if ~isfield(betaSummary, 'sessionRecords')
        error('betaSummary does not contain sessionRecords.');
    end

    if ~isfield(acrossOffsetSummary, 'fileInfo')
        error('acrossOffsetSummary does not contain fileInfo.');
    end

    if ~isfield(betaSummary.sessionRecords, 'fileName')
        error('betaSummary.sessionRecords does not contain fileName.');
    end

    if ~istable(acrossOffsetSummary.fileInfo) || ...
            ~ismember('fileName', acrossOffsetSummary.fileInfo.Properties.VariableNames)
        error('acrossOffsetSummary.fileInfo must be a table with a fileName column.');
    end

    betaFileNames = {betaSummary.sessionRecords.fileName}';
    kernelFileNames = acrossOffsetSummary.fileInfo.fileName;

    betaFileNames = normalizeFileNameList(betaFileNames, true);
    kernelFileNames = normalizeFileNameList(kernelFileNames, false);

    betaFileNames = unique(betaFileNames, 'stable');
    kernelFileNames = unique(kernelFileNames, 'stable');

    betaKeys = lower(betaFileNames);
    kernelKeys = lower(kernelFileNames);

    betaOnlyMask = ~ismember(betaKeys, kernelKeys);
    kernelOnlyMask = ~ismember(kernelKeys, betaKeys);

    betaOnly = betaFileNames(betaOnlyMask);
    kernelOnly = kernelFileNames(kernelOnlyMask);

    fprintf('Beta summary files:   %d\n', numel(betaFileNames));
    fprintf('Kernel summary files: %d\n', numel(kernelFileNames));

    if isempty(betaOnly) && isempty(kernelOnly)
        fprintf('\nBeta and kernel summaries are consistent: all files match.\n\n');
        return;
    end

    fprintf('\nMismatch detected.\n');

    if ~isempty(betaOnly)
        fprintf('\nFiles in IDR_BetaSummary but not IDR_KernelSummary:\n');
        printIndentedList(betaOnly);
    end

    if ~isempty(kernelOnly)
        fprintf('\nFiles in IDR_KernelSummary but not IDR_BetaSummary:\n');
        printIndentedList(kernelOnly);
    end

    fprintf('\n');
end


function fileNames = normalizeFileNameList(fileNames, isBetaSummary)
% Convert file names to comparable base names.

    if isstring(fileNames)
        fileNames = cellstr(fileNames);
    elseif ischar(fileNames)
        fileNames = {fileNames};
    elseif iscell(fileNames)
        % OK
    else
        error('File names must be a cell array, string array, or char array.');
    end

    fileNames = fileNames(:);

    for iFile = 1:numel(fileNames)
        thisName = fileNames{iFile};

        if isstring(thisName)
            thisName = char(thisName);
        end

        [~, stem, ext] = fileparts(thisName);

        if isempty(ext)
            normalizedName = stem;
        else
            normalizedName = [stem ext];
        end

        if isBetaSummary
            normalizedName = regexprep(normalizedName, ...
                '_scalarNoiseRegression(?=\.mat$|$)', '', ...
                'ignorecase');
        end

        % Normalize probe capitalization so Probe10 and probe10 compare equal
        % while preserving the original display form after matching.
        fileNames{iFile} = normalizedName;
    end
end


function printIndentedList(fileNames)

    for iFile = 1:numel(fileNames)
        fprintf('  %s\n', fileNames{iFile});
    end
end