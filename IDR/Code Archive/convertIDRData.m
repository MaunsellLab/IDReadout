function convertIDRData()
%   Convert Lablib .dat files to .mat format
%   convertIDRData() uses the default folder specified by domainFolder(mfilename('fullpath')).
%   convertIDRData(dataFolder) uses the specified folder.
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'DatFiles'));
if ~existed 
  fprintf('  convertIDRData -- failed to find data to convert in %s', dataFolder);
  return;
end
[names, paths] = getMatFileList(fullfile('Data', 'DatFiles'), "dat");
convertedFolder = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'Converted'));

% Convert any unconverted .dat files
numSkipped = 0;
for fi = 1:numel(paths)
  datPath = string(paths(fi));
  datName = string(names(fi));
  [~, baseName] = fileparts(datPath);
  tempInfoFileName = sprintf('%s/%s_fileInfo.mat', dataFolder, baseName);
  infoFileName = sprintf('%s/%s_fileInfo.mat', convertedFolder, baseName); % .mat headerfile name for .dat file
  outFileName = sprintf('%s/%s.mat', convertedFolder, baseName);           % .mat data file name for .dat file
  if isfile(infoFileName) && isfile(outFileName)                    % Skip if *_fileInfo.mat already exists
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
  fprintf('       converting %s\n', datName);
  for t = 1:nTrials
    trials{t} = readLLFile('t', t);
    if t == nTrials || mod(t, 500) == 0                            % Occasionally update the waitbar text
      fprintf('       reading %s: trial %d of %d (%.0f%%)\n', datName, t, nTrials, 100 * t / nTrials);
    end
  end

  % correct indexing errors that existed in earlier version of IDR
  trials = correctIndices(header, trials);

  % Save header and trials; use -v7.3 if these can get large
  fprintf('Saving %s\n', outFileName);
  save(outFileName, 'trials', 'header');
end
% if numSkipped > 0
%   fprintf('     convertIDRData: Skipped %d previously converted files\n', numSkipped);
% end
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
