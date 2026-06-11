function rootPath = folderPath()
% folderPath  Return the repository root path.
%
% Assumes this file is two folder levels below the repository root.

rootPath = ancestorFolder(fileparts(mfilename('fullpath')), 2);
end

%%-- unpack folder depth

function p = ancestorFolder(p, n)
for k = 1:n
    p = fileparts(p);
end
end
