function convertIDRData()
% convertIDRData  Convert Lablib .dat files to .mat format if needed.
%   convertIDRData() uses the default dataFolder.
%   convertIDRData(dataFolder) uses the specified folder.
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

  % specifier = folderPath() + "/DatFiles";               % default to a list from the Data folder
  % if nargin > 0
  %   specifier = specifier + '/' + singleFile;
  % end
  path = char(folderPath());
  convertedFolder = fullfile(path, 'Data', 'Converted');
  dataFolder   = fullfile(path, 'Data', 'DatFiles');
  [names, paths] = getMatFileList(fullfile('Data', 'DatFiles'), "dat");
  
  % Convert any unconverted .dat files
  numSkipped = 0;
  for fi = 1:numel(paths)
    datPath = string(paths(fi));
    datName = string(names(fi));
    [~, baseName] = fileparts(datPath);
    tempInfoFileName = sprintf('%s/%s_fileInfo.mat', dataFolder, baseName);
    infoFileName = sprintf('%s/%s_fileInfo.mat', convertedFolder, baseName); % Expected .mat headerfile name for this .dat file
    outFileName = sprintf('%s/%s.mat', convertedFolder, baseName);           % Expected .mat data file name for this .dat file
    if isfile(infoFileName) && isfile(outFileName)                    % Skip if *_fileInfo.mat already exists
      % fprintf('convertIDRData: Skipping %s (output exists)\n', datName);
      numSkipped = numSkipped + 1;
      continue;
    end
    if isfile(infoFileName) && ~isfile(outFileName)                   % Skip only *_fileInfo.mat exists, re-read
      delete(infoFileName);
    end
  
    fprintf('Reading %s\n', datName);
    header = readLLFile('i', datPath);
    movefile(tempInfoFileName, infoFileName);   % readLLFile leaves the info with the dat -- move it to the mat
    nTrials = header.numberOfTrials;
    trials  = cell(1, nTrials);
    fprintf('Converting %s\n', datName);
    for t = 1:nTrials
      trials{t} = readLLFile('t', t);   
      if t == nTrials || mod(t, 500) == 0                            % Occasionally update the waitbar text
        fprintf('Reading %s: Trial %d of %d (%.0f%%)\n', datName, t, nTrials, 100 * t / nTrials);
      end
    end

    % correct indexing errors that existed in earlier version of IDR
    trials = correctIndices(header, trials);

    % Save header and trials; use -v7.3 if these can get large
    fprintf('Saving %s', outFileName);
    save(outFileName, 'trials', 'header');
    % fprintf('Converted %s.dat → %s.mat & %s_fileInfo.mat\n', baseName, baseName, baseName);
  end
  if numSkipped > 0
    fprintf(' convertIDRData: Skipped %d previously converted files\n', numSkipped);
  end
end

% Several events were mislabeled as "change..." and "noChange..." when they
% were in fact "RF..." and "Opp...". These include "...Dots",
% "...PrefCohsPC", "...ProbeCohsPC", and "...TimesMS".  Specifically, the RF location
% was always assigned to be the "change" side and the oppposite location was alway 
% assigned to be the noChange side. Thus, whenever the change was on the opposite side,
% the assignments were backwards.  Here, we patch the change/noChange assignments 
% to correct this error by reassigning the values.

function trials = correctIndices(header, trials)

d = datetime(strtrim(string(header.date)), 'InputFormat', 'MMMM d, yyyy');
if d < datetime(2026, 2, 13) || d > datetime(2026, 3, 12)
    return;
else
  fprintf('Data from %s. Correcting trial-shift error.\n', header.date);
  nTrials  = numel(trials);
  for k = 1:nTrials
    tr = trials{k};
    if ~isfield(tr, 'trialEnd')
      continue;                       % skip if missing outcome
    end

    if k < nTrials
      patchTrial = trials{k+1};
    else
      patchTrial = tr;
    end
    % Trial shift without swapping
    tr.changeDots =           patchTrial.changeDots;
    tr.noChangeDots =         patchTrial.noChangeDots;
    tr.changePrefCohsPC =     patchTrial.changePrefCohsPC;
    tr.changeProbeCohsPC =    patchTrial.changeProbeCohsPC;
    tr.noChangePrefCohsPC =   patchTrial.noChangePrefCohsPC;
    tr.noChangeProbeCohsPC =  patchTrial.noChangeProbeCohsPC;
    tr.changeTimesMS =        patchTrial.changeTimesMS;
    tr.noChangeTimesMS =      patchTrial.noChangeTimesMS;
    trials{k} = tr;
  end
end
end
