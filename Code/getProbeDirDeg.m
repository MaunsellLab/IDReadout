function probeDirDeg = getProbeDirDeg(filePath)

S = load(filePath, 'header');

if ~isfield(S, 'header') || ~isfield(S.header, 'probeDirDeg')
    error('Missing header.probeDirDeg in %s', filePath);
end

probeDirDeg = S.header.probeDirDeg.data;

if ~isscalar(probeDirDeg) || ~isnumeric(probeDirDeg)
    error('Invalid probeDirDeg in %s', filePath);
end
end