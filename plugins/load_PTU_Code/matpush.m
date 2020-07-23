function matpush(mfile,varname,values,dim)
% Appends VALUES to variable named VARNAME within the matfile MFILE.
% If VARNAME or MFILE do not exist they are created.
% It is not possible to add extra dimensions.
%
% Jan Christoph Thiele, christoph.thiele@phys.uni-goettingen.de, 2020
%
narginchk(3,4);
mf = matfile(mfile,'Writable',true);

if ismember(varname,who(mf)) % var exists
    vsize = size(mf,varname);
    if prod(vsize)>0
        if nargin>3 && ~isempty(dim)
            if dim<1 && dim>numel(vsize)
                error('Dimensions %g are not supported.',dim);
            end
        else % First non empty dimension in the existing matfile
            dim = find(vsize~=1, 1 );
            if isempty(dim), dim = 1; end
        end
        
        if size(values,dim) > 0
            s = struct;
            s(1).type = '.';
            s(1).subs = varname;
            s(2).type = '()';
            s(2).subs = repmat({':'},numel(vsize),1);
            s(2).subs(vsize==0) = {1};
            s(2).subs(dim) = {vsize(dim)+1:vsize(dim)+size(values,dim)};
            
            subsasgn(mf,s,values); %#ok<SUBSASGN> %No output needed when wokring on matfiles.
        end
    else
        mf.(varname) = values;
    end
else
    mf.(varname) = values;
end


