
function updateIDS(varargin)
% updateIDS  Run the IDS analysis pipeline.
%
% convertIDSData is expected to create Data/FullSessions/*.mat. The remaining
% stages are supplied with this analysis package.

fprintf('>>> dailyUpdate start\n');
animal = 'All';
includeIDQ = parameterValue(varargin, 'IncludeIDQ', true);

fprintf('  >> convertIDSData start\n');
convertIDSData();
fprintf('  << convertIDSData complete\n');

fprintf('  >> makeIDSProbeSessions start\n');
makeIDSProbeSessions(varargin{:});
fprintf('  << makeIDSProbeSessions complete\n');

fprintf('  >> makeIDSSessionAnalyses start\n');
makeIDSSessionAnalyses(varargin{:});
fprintf('  << makeIDSSessionAnalyses complete\n');

fprintf('  >> plotIDSSummary start\n');
aSummaries = plotIDSSummary('Animal', 'All', 'IncludeIDQ', includeIDQ);
plotIDSGainAcrossOffsets(aSummaries);
mSummaries = plotIDSSummary('Animal', 'Meetz', 'IncludeIDQ', includeIDQ);
plotIDSGainAcrossOffsets(mSummaries);
nSummaries = plotIDSSummary('Animal', 'Neesha', 'IncludeIDQ', includeIDQ);
plotIDSGainAcrossOffsets(nSummaries);
fprintf('  << plotIDSSummary complete\n');

fprintf('<<< dailyUpdate complete\n');
end

function value = parameterValue(args, parameterName, defaultValue)
value = defaultValue;
for k = 1:2:numel(args)
  if (ischar(args{k}) || isstring(args{k})) && strcmpi(char(args{k}), parameterName)
    assert(k < numel(args), 'updateIDS:MissingParameterValue', ...
      'Missing value for parameter %s.', parameterName);
    value = args{k + 1};
    return
  end
end

end