function x = get(S,names,fallback)
%GET First available field, with Lablib .data wrappers removed.
if ischar(names)||isstring(names), names=cellstr(names); end
x=fallback;
for k=1:numel(names)
  if isstruct(S)&&isfield(S,names{k})
    x=mtexport.value(S.(names{k})); return
  end
end
end
