function [type, castedvar] =  intminclass(nvar, unsigned_flag)
% function [type, castedvar] =  minintclass(var)
%  Returns the most compact integer class which can accomodate var.
%  Second output contains the converted var.

if nargin<2 || isempty(unsigned_flag)
    if min(nvar(:))>=0
        unsigned_flag = true;
    else
        unsigned_flag = false;
    end
end

if unsigned_flag
    if all(nvar(:)<intmax('uint8')) || isempty(nvar)
        type = 'uint8';
    elseif all(nvar(:)<intmax('uint16'))
        type = 'uint16';
    elseif all(nvar(:)<intmax('uint32'))
        type = 'uint32';
    elseif all(nvar(:)<intmax('uint64'))
        type = 'uint64';
    else
        error('Variabel to large for uint64');
    end
else
    if all(nvar(:)<intmax('int8') & nvar(:)>intmin('int8')) || isempty(nvar)
        type = 'int8';
    elseif all(nvar(:)<intmax('int16') & nvar(:)>intmin('int16'))
        type = 'int16';
    elseif all(nvar(:)<intmax('int32') & nvar(:)>intmin('int32'))
        type = 'int32';
    elseif all(nvar(:)<intmax('int64') & nvar(:)>intmin('int64'))
        type = 'int64';
    else
        error('Variabel to large for int64');
    end
end

if nargout > 1
    castedvar = cast(nvar,type);
end