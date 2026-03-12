function convertIDRData(singleFile)
% convertIDRData  Convert Lablib .dat files to .mat format if needed.
%   convertIDRData() uses the default dataFolder.
%   convertIDRData(dataFolder) uses the specified folder.
%
%   For each  .dat file in dataFolder, this function checks for a
%   corresponding *_fileInfo.mat file. If it does not exist, the .dat
%   file is read with readLLFile and the trials/header are saved.

specifier = "Data";               % default to a list from the Data folder
if nargin > 0
  specifier = specifier + '/' + singleFile;
end
[names, paths] = getMatFileList(specifier, "dat");
  
  % Convert any unconverted .dat files
  for fi = 1:numel(paths)
    datPath = string(paths(fi));
    datName = string(names(fi));
    [filePath, baseName] = fileparts(datPath);
    infoFileName = sprintf('%s/%s_fileInfo.mat', filePath, baseName); % Expected .mat headerfile name for this .dat file
    outFileName = sprintf('%s/%s.mat', filePath, baseName);           % Expected .mat data file name for this .dat file
    if isfile(infoFileName) && isfile(outFileName)                    % Skip if *_fileInfo.mat already exists
      fprintf('convertIDRData: Skipping %s (output exists)\n', datName);
      continue;
    end
    if isfile(infoFileName) && ~isfile(outFileName)                   % Skip only *_fileInfo.mat exists, re-read
      delete(infoFileName);
    end
  
    % Only proceed if the .dat file itself is present (paranoid check)
    if ~isfile(datPath)
      fprintf('Warning: %s not found, skipping.\n', datName);
      continue;
    end
    waitTitle = sprintf('Reading %s', datName);
    f = waitbar(0, waitTitle);
    header = readLLFile('i', datPath);
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

    % correct indexing errors that existed in earlier version of IDR
    trials = correctIndices(header, trials);

    % Save header and trials; use -v7.3 if these can get large
    save(outFileName, 'trials', 'header');
    close(f);
    fprintf('Converted %s → %s & %s\n', datPath, infoFileName, outFileName);
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

versionStr = header.text.data{3};
version = str2double(regexp(versionStr,'\d+(\.\d+)?$','match','once'));
if (version >= 2.0)
  fprintf('IDR Data Version is %.1f. No index correction required\n', version);
  return;
else
  fprintf('IDR Data Version is %.1f. Correcting indexes\n', version);
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
    % If the change was on the RF side, frame shift without swapping
    if tr.trial.data.changeSide == 0
      tr.changeDots =           patchTrial.changeDots;
      tr.noChangeDots =         patchTrial.noChangeDots;
      tr.changePrefCohsPC =     patchTrial.changePrefCohsPC;
      tr.changeProbeCohsPC =    patchTrial.changeProbeCohsPC;
      tr.noChangePrefCohsPC =   patchTrial.noChangePrefCohsPC;
      tr.noChangeProbeCohsPC =  patchTrial.noChangeProbeCohsPC;
      tr.changeTimesMS =        patchTrial.changeTimesMS;
      tr.noChangeTimesMS =      patchTrial.noChangeTimesMS;
      
    % Otherwise change on Opp side and we need to swap the assignments
    else
      tr.changeDots =           patchTrial.noChangeDots;
      tr.noChangeDots =         patchTrial.changeDots;
      tr.changePrefCohsPC =     patchTrial.noChangePrefCohsPC;
      tr.changeProbeCohsPC =    patchTrial.noChangeProbeCohsPC;
      tr.noChangePrefCohsPC =   patchTrial.changePrefCohsPC;
      tr.noChangeProbeCohsPC =  patchTrial.changeProbeCohsPC;
      tr.changeTimesMS =        patchTrial.noChangeTimesMS;
      tr.noChangeTimesMS =      patchTrial.changeTimesMS;
    end

    trials{k} = tr;
  end
end
end

% function trials = correctIndices(header, trials)
% 
% versionStr = header.text.data{3};
% version = str2double(regexp(versionStr,'\d+(\.\d+)?$','match','once'));
% if (version >= 2.0)
%   fprintf('IDR Data Version is %.1f. No index correction required\n', version);
%   return;
% else
%   fprintf('IDR Data Version is %.1f. Correcting indexes\n', version);
%   nTrials  = numel(trials);
%   for k = 1:nTrials
%     tr = trials{k};
%     if ~isfield(tr, 'trialEnd')
%       continue;                       % skip if missing outcome
%     end
%     % If the change was on the RF side, everything is correct
%     changeSide = tr.trial.data.changeSide;
%     if changeSide == 0               % change on Opposite side -- misassigned
%       continue;
%     end
% 
%     % Otherwise we need to swap the assignments
%     thisTrial = tr;
%     tr.changeDots =           thisTrial.noChangeDots;
%     tr.noChangeDots =         thisTrial.changeDots;
%     tr.changePrefCohsPC =     thisTrial.noChangePrefCohsPC;
%     tr.changeProbeCohsPC =    thisTrial.noChangeProbeCohsPC;
%     tr.noChangePrefCohsPC =   thisTrial.changePrefCohsPC;
%     tr.noChangeProbeCohsPC =  thisTrial.changeProbeCohsPC;
%     tr.changeTimesMS =        thisTrial.noChangeTimesMS;
%     tr.noChangeTimesMS =      thisTrial.changeTimesMS;
% 
%     trials{k} = tr;
%   end
% end
% end