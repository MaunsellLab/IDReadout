%% validFolder()
function [fullPath, existed] = validFolder(fullPath)
% check whether a folder exists and create it otherwise

  if ~exist(fullPath, 'dir')
    mkdir(fullPath);
    existed = false;
  else
    existed = true;
  end
end
