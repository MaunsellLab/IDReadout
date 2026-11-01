function updateIDL(varargin)
% updateIDL  Update sessions, then analyse each IDL allocation style separately.
% updateIDL('Animal','Meetz') processes every available style independently.
% updateIDL('Animal','Meetz','TwoAllocationMode',true) selects the new style.
p=inputParser;
p.KeepUnmatched=true;
p.PartialMatching=false;
p.addParameter('Animal','All',@(x) ischar(x)||isstring(x));
p.addParameter('TwoAllocationMode',[],@(x) isempty(x)||(islogical(x)&&isscalar(x)));
p.addParameter('Replace',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('ReplaceAnalysis',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('NBootstrap',1000,@(x) isnumeric(x)&&isscalar(x)&&x>=0&&x==round(x));
p.addParameter('IncludeNoiseFreeSessions',false,@(x) islogical(x)&&isscalar(x));
p.addParameter('Visible','on',@(x) any(strcmpi(string(x),["on" "off"])));
p.addParameter('Verbose',true,@(x) islogical(x)&&isscalar(x));
p.parse(varargin{:});
o=p.Results;
names=fieldnames(p.Unmatched); selection=cell(1,2*numel(names));
for k=1:numel(names)
  selection{2*k-1}=names{k}; selection{2*k}=p.Unmatched.(names{k});
end
fprintf('>>> updateIDL start\n');
convertIDLData();
makeIDLSessionAnalyses('Animal',o.Animal,'TwoAllocationMode',o.TwoAllocationMode, ...
  'Replace',o.ReplaceAnalysis,'Verbose',o.Verbose,selection{:});
makeIDLSessionSummaryPlots('Animal',o.Animal,'TwoAllocationMode',o.TwoAllocationMode, ...
  'Replace',o.Replace,'Verbose',o.Verbose,selection{:});
root=domainFolder(mfilename('fullpath'));
folder=fullfile(root,'Data','SessionAnalysis');
modes=[false true];
if ~isempty(o.TwoAllocationMode), modes=o.TwoAllocationMode; end
for mode=modes
  [files,info]=selectAnalysisFiles(folder,'Animal',o.Animal, ...
    'TwoAllocationMode',mode,selection{:});
  if isempty(files), continue, end
  fprintf('  IDL allocation style: %s\n',styleName(mode));
  hasNoise=false;
  for s=1:numel(files)
    S=load(files{s},'sessionAnalysis'); A=S.sessionAnalysis;
    hasNoise=hasNoise || any(abs(double(A.changeNoiseByFrameTrial))>1e-8,'all') || ...
      any(abs(double(A.noChangeNoiseByFrameTrial))>1e-8,'all');
  end
  if hasNoise || o.IncludeNoiseFreeSessions
    [model,~]=fitIDLModelRecoveryAcrossSessions('Animal',o.Animal, ...
      'TwoAllocationMode',mode,'IncludeNoiseFreeSessions',o.IncludeNoiseFreeSessions, ...
      'Visible',o.Visible,selection{:});
    plotIDLAcrossSessionBehavior('Animal',o.Animal,'TwoAllocationMode',mode, ...
      'IncludeNoiseFreeSessions',o.IncludeNoiseFreeSessions,'ModelResult',model, ...
      'Visible',o.Visible,selection{:});
  else
    fprintf('  Noise-model comparison skipped: no recorded step noise.\n');
  end
  animals=unique(string(info.animal),'stable');
  for animal=animals(:)'
    n=sum(strcmpi(string(info.animal),animal));
    if n<3
      fprintf('  p-norm profile skipped for %s: %d sessions; at least 3 required.\n',animal,n);
      continue
    end
    pilotIDLPNormProfile('Animal',char(animal),'TwoAllocationMode',mode, ...
      'NBootstrap',o.NBootstrap,'Verbose',o.Verbose,'Visible',o.Visible,selection{:});
  end
end
fprintf('<<< updateIDL complete\n');
end

function name=styleName(two)
name='FourAllocations';
if two, name='TwoAllocations'; end
end
