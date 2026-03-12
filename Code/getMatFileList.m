function [names, paths] = getMatFileList(target, extension)
% getMatFileList Return file names for a given exteions and full paths from standard folders
%
% Usage:
%   [names, paths] = getMatFileList('Data', 'dat')
%   [names, paths] = getMatFileList('Kernels')
%   [names, paths] = getMatFileList('Data/session01')
%
% Input:
%   target    string or char
%
% Output:
%   names     cell array of file names with specified extention (defaults to .mat)
%   paths     cell array of full file paths for those files

  if nargin < 2
    wildcardExtension = sprintf('*.mat');
  else
    wildcardExtension = sprintf('*.%s', extension);
  end

    validSubfolders = {'Data','Kernels','Models','Plots','Summaries'};
    rootFolder = folderPath();

    if isstring(target)
        target = char(target);
    end

    target = strtrim(target);
    target = strrep(target,'\','/');   % normalize separators

    % ---------- CASE 1: entire subfolder ----------
    if any(strcmp(target, validSubfolders))

        subfolder = target;
        folder = fullfile(rootFolder, subfolder);

        if ~isfolder(folder)
            error('Folder does not exist: %s', folder);
        end

        S = dir(fullfile(folder, wildcardExtension));

        names = {S.name}';
        paths = fullfile(folder, names);

        % optional sorting
        [names, idx] = sort(names);
        paths = paths(idx);

        return
    end

    % ---------- CASE 2: single file ----------
    parts = split(string(target),'/');
    parts = parts(parts ~= "");

    if numel(parts) ~= 2
        error('Input must be "Subfolder" or "Subfolder/filename".');
    end

    subfolder = char(parts(1));
    fileStem  = char(parts(2));

    if ~any(strcmp(subfolder, validSubfolders))
        error('Invalid subfolder: %s', subfolder);
    end

    if ~endsWith(fileStem, extension)
        % fileName = [fileStem extension];
        fileName = sprintf('%s.%s', fileStem , extension);
    else
        fileName = fileStem;
    end

    names = {fileName};
    paths = {fullfile(rootFolder, subfolder, fileName)};
end