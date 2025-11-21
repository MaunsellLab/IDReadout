function convertIDPhysTraining()
    % convert dat files to mat format
    % Requires readLLFile
    
    dataFolder = '/Users/maunsell/Desktop/IDR/Data';
    cd(dataFolder);

    datFileList = dir('*.dat');                                             % All the .dat files
    matFileList = dir('*_fileInfo.mat');                                    % All the converted .mat files
                                                                            % Consider only the header files just for making
                                                                            % judgments about whether the .dat file needs to
                                                                            % be converted
    
    for fi = 1:length(datFileList)
        if fi - length(matFileList) > 0
            if isfile(datFileList(fi).name)
                   f = waitbar(0, 'Reading the Headers');
                   
                   header = readLLFile('i', datFileList(fi).name);
                   trials = cell(1,header.numberOfTrials);
                   waitbar(0.2, f, ['Reading the Trials from ', datFileList(fi).name])
                   pause(2);
                   nTrials = header.numberOfTrials;
                   for i = 1:nTrials
                       trials{i} = readLLFile('t', i);
                       waitbar((0.2 + i / nTrials * 0.69), f, ...
                           sprintf(['Reading the Trials from ', datFileList(fi).name, '\nProgress %d %%'], ...
                           floor((0.2 + i / nTrials * 0.69)*100)))
                   end
                   waitbar(0.99, f, 'Saving the .mat Files')
                   matFileName = strrep(datFileList(fi).name, '.dat', '.mat');
                   save(matFileName, 'trials', 'header')
                   close(f)
            end
        end
    end
end