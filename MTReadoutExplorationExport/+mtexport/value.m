function x = value(x)
%VALUE Unwrap Lablib event/header data without modifying the source.
while isstruct(x)&&isscalar(x)&&isfield(x,'data'), x=x.data; end
end
