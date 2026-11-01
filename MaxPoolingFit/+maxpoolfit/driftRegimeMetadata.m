function [sourceDesign,trialGroup] = driftRegimeMetadata(cache,varargin)
%DRIFTREGIMEMETADATA Classify trials by drift-direction randomization.
%
% IDQ sources are classified from sessionHeader.randomBaseDir when that
% field is available. Whether cached drift directions depart from the fixed
% 30/150/270-degree triplet is an independent audit and is the fallback for
% historical or unavailable headers. IDR and IDS retain their domain labels.

p=inputParser;
p.addParameter('ReadHeaders',true,@(x)islogical(x)&&isscalar(x));
p.addParameter('HeaderField',"randomBaseDir",@(x)isstring(x)||ischar(x));
p.addParameter('FixedBaseDirectionsDeg',[30 150 270], ...
  @(x)isnumeric(x)&&isvector(x)&&~isempty(x)&&all(isfinite(x)));
p.addParameter('DirectionToleranceDeg',0.1, ...
  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
p.parse(varargin{:});
maxpoolfit.validateCache(cache);
T=cache.trials;
n=height(T);

if ismember('sourceFile',T.Properties.VariableNames)
  sourceFile=string(T.sourceFile);
else
  sourceFile=strings(n,1);
end
key=string(T.domain)+"|"+string(T.animal)+"|"+string(T.sessionID)+"|"+sourceFile;
[~,~,sourceIndex]=unique(key);
nSources=max(sourceIndex);
firstRow=accumarray(sourceIndex,(1:n)',[nSources 1],@min);

domain=string(T.domain(firstRow));
animal=string(T.animal(firstRow));
sessionID=string(T.sessionID(firstRow));
sourceFile=sourceFile(firstRow);
nCachedTrials=accumarray(sourceIndex,1,[nSources 1]);
nUniqueDriftDirections=zeros(nSources,1);
randomBaseDir=nan(nSources,1);
realizedRandomBaseDir=false(nSources,1);
headerAvailable=false(nSources,1);
classificationSource=strings(nSources,1);
analysisGroup=strings(nSources,1);
auditMismatch=false(nSources,1);
headerField=char(string(p.Results.HeaderField));

for k=1:nSources
  idx=find(sourceIndex==k);
  d=mod(double(T.driftDirectionDeg(idx)),360);
  d=round(d.*1e6)./1e6;
  nUniqueDriftDirections(k)=numel(unique(d(isfinite(d))));

  switch domain(k)
    case "IDR"
      analysisGroup(k)="IDR";
      classificationSource(k)="domain";
    case "IDQ"
      fixedDirs=mod(double(p.Results.FixedBaseDirectionsDeg(:)).',360);
      circularDistance=abs(mod(d-fixedDirs+180,360)-180);
      realizedRandomBaseDir(k)=any(min(circularDistance,[],2)> ...
        p.Results.DirectionToleranceDeg);
      [headerAvailable(k),randomBaseDir(k)]=readHeaderFlag( ...
        sourceFile(k),headerField,p.Results.ReadHeaders);
      if headerAvailable(k)
        isRandom=randomBaseDir(k)~=0;
        classificationSource(k)="sessionHeader."+string(headerField);
        auditMismatch(k)=isRandom~=realizedRandomBaseDir(k);
      else
        isRandom=realizedRandomBaseDir(k);
        classificationSource(k)="observed driftDirectionDeg";
      end
      if isRandom
        analysisGroup(k)="IDQ-randomized";
      else
        analysisGroup(k)="IDQ-fixed";
      end
    case "IDS"
      analysisGroup(k)="IDS";
      classificationSource(k)="domain";
    otherwise
      error('maxpoolfit:driftRegimeMetadata:UnknownDomain', ...
        'Cannot classify domain %s.',domain(k));
  end
end

sourceDesign=table(sourceFile,sessionID,domain,animal,nCachedTrials, ...
  randomBaseDir,realizedRandomBaseDir,nUniqueDriftDirections, ...
  headerAvailable,classificationSource,auditMismatch,analysisGroup);
trialGroup=analysisGroup(sourceIndex);

if any(auditMismatch)
  warning('maxpoolfit:driftRegimeMetadata:HeaderDirectionMismatch', ...
    ['%d IDQ source(s) have sessionHeader.%s values inconsistent with ' ...
     'the cached drift directions relative to the fixed base-direction ' ...
     'triplet. The header value ' ...
     'was retained; inspect sourceDesign(sourceDesign.auditMismatch,:).'], ...
    sum(auditMismatch),headerField);
end
end

function [available,value]=readHeaderFlag(sourceFile,fieldName,readHeaders)
available=false; value=NaN;
if ~readHeaders || strlength(sourceFile)==0 || ~isfile(sourceFile)
  return
end
S=load(char(sourceFile),'sessionHeader');
if ~isfield(S,'sessionHeader') || ~isstruct(S.sessionHeader) || ...
    ~isfield(S.sessionHeader,fieldName)
  return
end
v=unwrapScalar(S.sessionHeader.(fieldName));
if islogical(v) || isnumeric(v)
  if isscalar(v) && isfinite(double(v))
    value=double(logical(v)); available=true;
  end
elseif ischar(v) || (isstring(v)&&isscalar(v))
  s=lower(strtrim(string(v)));
  if ismember(s,["true" "yes" "on" "1"])
    value=1; available=true;
  elseif ismember(s,["false" "no" "off" "0"])
    value=0; available=true;
  end
end
end

function v=unwrapScalar(v)
while iscell(v)&&isscalar(v), v=v{1}; end
if isstruct(v)&&isscalar(v)&&isfield(v,'data'), v=unwrapScalar(v.data); end
end
