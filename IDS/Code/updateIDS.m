function updateIDS(varargin)
% updateIDS  Run the IDS analysis pipeline.
%
% convertIDSData is expected to create Data/FullSessions/*.mat. The remaining
% stages are supplied with this analysis package.

convertIDSData();
makeIDSProbeSessions(varargin{:});
makeIDSSessionAnalyses(varargin{:});
plotIDSSummary();
end
