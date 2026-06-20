function values = getTrialField(trials, varargin)
% getTrialField
%
% Extract values from a cell array of trial structs.
%
% Examples:
%   eotCode = getTrialField(trials, 'extendedEOT', 'data');
%   sideIdx = getTrialField(trials, 'sideIndex');
%   x       = getTrialField(trials, {'extendedEOT', 'data'});
%
% Missing fields intentionally crash.

if isscalar(varargin) && iscell(varargin{1})
  fieldPath = varargin{1};
else
  fieldPath = varargin;
end

if ~iscell(trials)
  trials = num2cell(trials);
end

nTrials = numel(trials);
values = cell(nTrials, 1);

for iTrial = 1:nTrials
  x = trials{iTrial};

  for iField = 1:numel(fieldPath)
    x = x.(fieldPath{iField});
  end

  values{iTrial} = x;
end

% Convert to numeric/logical array when possible.
if all(cellfun(@(x) isnumeric(x) || islogical(x), values)) && ...
    all(cellfun(@isscalar, values))
  values = cell2mat(values);
end

end