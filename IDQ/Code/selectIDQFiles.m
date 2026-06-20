function files = selectIDQFiles(folderPath, varargin)
% selectIDQFiles
%
% Centralized selector for IDQ converted session files.
%
% First-pass behavior:
%   - return *.mat files
%   - exclude *_fileInfo.mat files
%
% Varargin is reserved for future file/session inclusion criteria.

% Keep varargin accepted now so analysis functions can define selectArgs at
% the top before all selection criteria are implemented.
unusedArgs = varargin; %#ok<NASGU>

d = dir(fullfile(folderPath, '*.mat'));

if isempty(d)
  files = struct('name', {}, 'path', {});
  return
end

names = {d.name};

exclude = endsWith(names, '_fileInfo.mat');
d = d(~exclude);

files = struct('name', {}, 'path', {});
for iFile = 1:numel(d)
  files(iFile).name = d(iFile).name;
  files(iFile).path = fullfile(d(iFile).folder, d(iFile).name);
end

end