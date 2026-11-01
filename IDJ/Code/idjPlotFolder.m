function folder = idjPlotFolder(varargin)
% idjPlotFolder  Return a folder below IDJ/Plots (source lives in IDJ/Code).
codeFolder = fileparts(mfilename('fullpath'));
idjRoot = fileparts(codeFolder);
folder = validFolder(fullfile(idjRoot, 'Plots', varargin{:}));
end
