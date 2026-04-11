function fileList = findSessionsByProbe(probeDirDeg, dataDir)

tag = probeDirTag(probeDirDeg);
probeDir = fullfile(dataDir, tag);

if ~exist(probeDir,'dir')
    fileList = {};
    return;
end

d = dir(fullfile(probeDir, '*.mat'));
fileList = cell(numel(d),1);

for i = 1:numel(d)
    filePath = fullfile(d(i).folder, d(i).name);

    p = getProbeDirDeg(filePath);
    if p ~= probeDirDeg
        error('Probe mismatch: file %s is %d, expected %d', ...
            d(i).name, p, probeDirDeg);
    end

    fileList{i} = filePath;
end
end