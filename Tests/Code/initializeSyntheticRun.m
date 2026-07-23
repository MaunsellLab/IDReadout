function manifest = initializeSyntheticRun(context,varargin)
% initializeSyntheticRun  Create one isolated, manifest-backed simulation.
%
% manifest = initializeSyntheticRun(context, ...
%   'ModelName',"MTMax", ...
%   'ModelParameters',parameters, ...
%   'SourceFiles',sourcePaths, ...
%   'CodeVersion',gitCommit)
%
% The run must not already contain files. Experimental SourceFiles are read
% only to record provenance; subsequent synthetic analyses are restricted
% to the run tree by their context.

context = validateAnalysisContext(context);
if context.Mode ~= "synthetic"
  error("IDReadout:SyntheticContextRequired", ...
    "initializeSyntheticRun requires a synthetic context.");
end

p = inputParser;
addParameter(p,'ModelName',"",@(x) ischar(x)||isstring(x));
addParameter(p,'ModelParameters',struct(),@(x) isstruct(x));
addParameter(p,'SourceFiles',strings(0,1),@(x) ischar(x)||isstring(x)||iscellstr(x));
addParameter(p,'CodeVersion',"",@(x) ischar(x)||isstring(x));
parse(p,varargin{:});
opts = p.Results;

runRoot = string(domainFolder(mfilename('fullpath'),'Tests',context));
if isfolder(runRoot)
  contents = dir(runRoot);
  contents = contents(~ismember({contents.name},{'.','..'}));
  if ~isempty(contents)
    error("IDReadout:SyntheticRunAlreadyExists", ...
      "Synthetic run is not empty and will not be reused: %s",runRoot);
  end
end

folders = string({ ...
  analysisPath(context,'IDQ','Data'); ...
  analysisPath(context,'IDQ','Data','AcrossSessionSummaries'); ...
  analysisPath(context,'IDQ','Plots'); ...
  analysisPath(context,'IDQ','Plots','AcrossSessionSummaries'); ...
  analysisPath(context,'IDR','Data'); ...
  analysisPath(context,'IDR','Data','AcrossSessionSummaries'); ...
  analysisPath(context,'IDR','Plots'); ...
  analysisPath(context,'IDR','Plots','AcrossSessionSummaries'); ...
  analysisPath(context,'Common Code','Data'); ...
  analysisPath(context,'Common Code','Data','AcrossSessionSummaries'); ...
  analysisPath(context,'Common Code','Plots'); ...
  analysisPath(context,'Common Code','Plots','AcrossSessionSummaries')});
for k = 1:numel(folders)
  if ~isfolder(folders(k))
    [ok,message] = mkdir(folders(k));
    if ~ok
      error("IDReadout:SyntheticRunCreateFailed", ...
        "Could not create %s: %s",folders(k),message);
    end
  end
end

sourceFiles = string(opts.SourceFiles(:));
sourceRecords = repmat(struct('Path',"",'Bytes',0,'Modified',"", ...
  'SHA256',""),numel(sourceFiles),1);
for k = 1:numel(sourceFiles)
  if ~isfile(sourceFiles(k))
    error("IDReadout:MissingSimulationSource", ...
      "Simulation source file does not exist: %s",sourceFiles(k));
  end
  info = dir(sourceFiles(k));
  sourceRecords(k).Path = canonicalPath(sourceFiles(k));
  sourceRecords(k).Bytes = info.bytes;
  sourceRecords(k).Modified = string(datetime(info.datenum, ...
    'ConvertFrom','datenum','TimeZone','local'));
  sourceRecords(k).SHA256 = fileSHA256(sourceFiles(k));
end

manifest = struct();
manifest.ManifestSchemaVersion = 1;
manifest.DataOrigin = "synthetic";
manifest.RunID = context.RunID;
manifest.CreatedAt = string(datetime('now','TimeZone','local'));
manifest.RandomSeed = context.Seed;
manifest.MaxRunBytes = context.MaxRunBytes;
manifest.SaveLatents = context.SaveLatents;
manifest.ModelName = string(opts.ModelName);
manifest.ModelParameters = opts.ModelParameters;
manifest.CodeVersion = string(opts.CodeVersion);
manifest.MATLABVersion = string(version);
manifest.ProjectRootAtCreation = context.ProjectRoot;
manifest.RunRootAtCreation = context.RunRoot;
manifest.SourceFiles = sourceRecords;
manifest.Status = "initialized";

manifestPath = analysisPath(context,'Tests','manifest.mat');
save(manifestPath,'manifest','-v7');
assertSimulationRunCapacity(context,0);
end

function digest = fileSHA256(path)
md = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(path,'rb');
if fid<0
  error("IDReadout:SourceReadFailed","Could not read source file: %s",path);
end
cleanup = onCleanup(@() fclose(fid));
while ~feof(fid)
  bytes = fread(fid,1024*1024,'*uint8');
  if ~isempty(bytes)
    md.update(typecast(bytes,'int8'));
  end
end
raw = typecast(md.digest(),'uint8');
digest = lower(string(reshape(dec2hex(raw,2).',1,[])));
clear cleanup
end

function path = canonicalPath(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end
