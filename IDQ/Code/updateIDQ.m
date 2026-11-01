function updateIDQ(varargin)
% updateIDQ

%% convertIDQData
% converts dat to mat

convertIDQData(varargin{:});

%% preprocessIDQSession 
% creates sessionHeader

preprocessIDQSession();

%% makeIDQAcrossSessionSummary 
% creates session summaries

plotIDQSummary();

