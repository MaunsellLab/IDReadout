function outStr = underscoreToDash(inStr)
% Replace underscores with dashes
if ischar(inStr)
  inStr = string(inStr);
end
outStr = replace(inStr, "_", "-");
end