function out = plain(in)
%PLAIN Serialize MATLAB analysis metadata without MCOS tables/strings.
if istable(in)
  out=struct('matlab_type','table','row_count',height(in));
  for k=1:width(in)
    key=in.Properties.VariableNames{k}; out.(key)=mtexport.plain(in.(key));
  end
elseif isstring(in)
  out=cellstr(in);
elseif isdatetime(in)||isduration(in)||iscategorical(in)
  out=cellstr(string(in));
elseif isa(in,'function_handle')
  out=func2str(in);
elseif isstruct(in)
  out=in;
  for j=1:numel(in)
    names=fieldnames(in(j));
    for k=1:numel(names), out(j).(names{k})=mtexport.plain(in(j).(names{k})); end
  end
elseif iscell(in)
  out=cellfun(@mtexport.plain,in,'UniformOutput',false);
elseif isnumeric(in)||islogical(in)||ischar(in)
  out=in;
else
  out=struct('matlab_type',class(in),'export_note','Unsupported object omitted');
end
end
