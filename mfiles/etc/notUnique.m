function NU = notUnique(A)
%%NOTUNIQUE -- get items that are not unique in A
%
% USAGE NU = notUnique(A)
%   A is array of class char, scalar or strings (cell array of strings)
%
% TO 140227 Deyang

uA = unique(A);

if numel(uA)==numel(A)
    NU = [];
else
    NU = false(1,numel(uA));
    switch class(A)
        case {'char','scalar'}
            for iu=1:numel(uA)
                if sum(uA(iu)==A)>1
                    NU(iu)=true;
                end
            end
        case 'cell'
            for iu=1:numel(uA)
                if numel(strmatchi(uA{iu},A))>1
                    NU(iu)=true;
                end
            end
        otherwise
            error('Argument should be of class char, scalar or cell array of strings, not <<%s>>!',class(A));
    end
    NU = uA(NU);
end