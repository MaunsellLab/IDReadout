function folder = idlPlotFolder(varargin)
% idlPlotFolder  Return an IDL plot folder relative to IDL/Code.
%
% IDL analysis source files live in IDL/Code. All graphical output is
% rooted at the sibling IDL/Plots folder.

codeFolder = fileparts(mfilename('fullpath'));
idlRoot = fileparts(codeFolder);
folder = validFolder(fullfile(idlRoot, 'Plots', varargin{:}));
end
