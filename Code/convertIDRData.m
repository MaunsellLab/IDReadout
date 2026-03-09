function convertIDRData(dataFolder)
% convertIDRData  Convert Lablib .dat files to .mat format if needed.
%   convertIDRData() uses the default dataFolder.
%   convertIDRData(dataFolder) uses the specified folder.
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

  if nargin < 1 || isempty(dataFolder)
    dataFolder = folderPath();
  end
  dataFolder = dataFolder + "/Data";
  
  % Remember original folder and restore it on exit
  originalFolder = pwd;
  cleanupObj = onCleanup(@() cd(originalFolder));
  cd(dataFolder);
  
  % List of all .dat files
  datFileList = dir('*.dat');
  if isempty(datFileList)
    fprintf('No .dat files found in %s\n', dataFolder);
    return;
  end
  
  % Convert any unconverted .dat files
  for fi = 1:numel(datFileList)
    datName = datFileList(fi).name;
    [~, baseName] = fileparts(datName);
    infoFileName = sprintf('%s_fileInfo.mat', baseName);             % Expected .mat headerfile name for this .dat file
    outFileName = sprintf('%s.mat', baseName);                       % Expected .mat data file name for this .dat file
    if isfile(infoFileName) && isfile(outFileName)                   % Skip if *_fileInfo.mat already exists
      fprintf('convertIDRData: Skipping %s (output exists)\n', datName, infoFileName, outFileName);
      continue;
    end
    if isfile(infoFileName) && ~isfile(outFileName)                   % Skip only *_fileInfo.mat exists, re-read
      delete(infoFileName);
    end
  
    % Only proceed if the .dat file itself is present (paranoid check)
    if ~isfile(datName)
      fprintf('Warning: %s not found, skipping.\n', datName);
      continue;
    end
    waitTitle = sprintf('Reading %s', datName);
    f = waitbar(0, waitTitle);
    header = readLLFile('i', datName);
    nTrials = header.numberOfTrials;
    trials  = cell(1, nTrials);
    for t = 1:nTrials
      trials{t} = readLLFile('t', t);   
      if t == nTrials || mod(t, 10) == 0                            % Occasionally update the waitbar text
        waitbar(t / nTrials, f, sprintf('Reading %s\nTrial %d of %d (%.0f%%)', ...
          datName, t, nTrials, 100 * t / nTrials));
      end
    end

    waitbar(1.0, f, sprintf('Saving %s', outFileName));
    % Save header and trials; use -v7.3 if these can get large
    save(outFileName, 'trials', 'header');
    close(f);
    fprintf('Converted %s → %s & %s\n', datName, infoFileName, outFileName);
  end
end
