function id = sessionID(name)
%SESSIONID Normalize equivalent derived/parent source-session identifiers.
[~,id]=fileparts(char(name));
id=regexprep(id,'_sessionAnalysis$','','ignorecase');
id=regexprep(id,'^ProbeSession','','ignorecase');
id=regexprep(id,'_Probe[-+]?\d+(?:\.\d+)?$','','ignorecase');
assert(~isempty(id),'Missing source-session identity.');
end
