function convertIDJData()
% convertIDLData  Convert IDL Lablib files into FullSession MAT files.
%
% Sources are in Data/DatFiles and outputs are in Data/FullSessions. The
% standalone makeIDLSessionHeader function supplies the analysis metadata.

[dataFolder, existed] = validFolder(fullfile(domainFolder(mfilename('fullpath')), 'Data', 'DatFiles'));
if ~existed
  fprintf('  convertIDJData -- failed to find data to convert in %s', dataFolder);
  return
end

S = dir(fullfile(dataFolder, '*.dat'));
names = {S.name}';
paths = fullfile(dataFolder, names);
[names, idx] = sort(names);
paths = paths(idx);
convertedFolder = validFolder(fullfile( ...
  domainFolder(mfilename('fullpath')), 'Data', 'FullSessions'));

for fileIndex = 1:numel(paths)
  datPath = string(paths(fileIndex));
  datName = string(names(fileIndex));
  [~, baseName] = fileparts(datPath);
  tempInfoPath = fullfile(dataFolder, baseName + '_fileInfo.mat');
  infoPath = fullfile(convertedFolder, baseName + '_fileInfo.mat');
  outputPath = fullfile(convertedFolder, baseName + '.mat');
  if isfile(infoPath) && isfile(outputPath)
    continue
  end
  if isfile(infoPath) && ~isfile(outputPath)
    delete(infoPath);
  end

  fprintf('      Reading %s\n', datName);
  header = readLLFile('i', datPath);
  movefile(tempInfoPath, infoPath);
  nTrials = header.numberOfTrials;
  trials = cell(1, nTrials);
  fprintf('       converting %s (%d trials)\n', datName, nTrials);
  for trialIndex = 1:nTrials
    trials{trialIndex} = readLLFile('t', trialIndex);
    if trialIndex == nTrials || mod(trialIndex, 500) == 0
      fprintf('       reading %s: trial %d of %d (%.0f%%)\n', ...
        datName, trialIndex, nTrials, 100 * trialIndex / nTrials);
    end
  end

  sessionHeader = makeIDJSessionHeader(header);
  save(outputPath, 'trials', 'header', 'sessionHeader');
  fprintf('      Saved %s\n', outputPath);
end
end
