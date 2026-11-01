function key = id(path)
%ID Portable cache filename; full fingerprints are separately checked.
b=double(unicode2native(char(path),'UTF-8')); a=5381; c=2166136261;
for k=1:numel(b)
  a=mod(a*33+b(k),4294967291); c=mod(c*65599+b(k),4294967279);
end
key=sprintf('%08x%08x',uint32(a),uint32(c));
end
