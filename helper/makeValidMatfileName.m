function name = makeValidMatfileName(name)
%makeValidMatfileName Sanitizes a char to make sure it can be used as matfile
%name. Unlike v7 matfiles v7.3 matfiles only accept ASCII filenames. This is not
%documented.
% Invalid characters are replaced with an urlencoded version. The characters
% \*/: are not escaped.
% If name is path this is applied to the whole path which can leed to a
% different output folder.

acceptedChar = '!#$%&''()+,-./0123456789:;=@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{}~';
valid = ismember(name, acceptedChar);
if ~all(valid)
    name = num2cell(name);
    name(~valid) = cellfun(@urlencode,name(~valid),'UniformOutput',false);
    name = cell2mat(name);
end
end

