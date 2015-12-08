function ltype=mf_linetype(it,ltype)
%MF_LINETYPE yields next in ltype series and restarts after last line type has been used
%
% Example:
%    ltype=mf_linetype(it[,ltypelist)
%    ltype=mf_linetype(it,{'-'  '--'  ':'  '-.'  'none'});
%    if nargin==1, ltype={'-' ';' };
%
% See alo: mf_color mf_marker
% 
% TO 120101 120531

linetypelist      ={'-'  '--'  ':'  '-.'};

if nargin==1
    NL=length(linetypelist);
    ic=rem(it,NL); if ic==0, ic=NL; end
    ltype=linetypelist{ic};
else 
    if ischar(ltype)
        if ~strmatchi('-',ltype)
            ltype='-';
        end
    elseif isnumeric(ltype),
        NL=ltype(1);
        if NL<=0,
            ltype={'-'};
        else
            ltype=linetypelist(1:min(NL,length(linetypelist)));
        end
    elseif iscell
        for ic=length(ltype:-1:1)
            if ~strmatchi(ltype{ic},linetypelist),
                ltype=ltype(1:ic-1);
            end
        end
        if isempty(ltype)
            ltype={'-'};
        end
    else
        ltype={'-'};
    end
    NL=length(ltype);
    ic=rem(it,NL); if ic==0, ic=NL; end
    ltype=ltype{ic};
end
            
% Matlab's line types list    
% '+' Plus sign
% 'o' Circle
% '*' Asterisk
% '.' Point
% 'x' Cross
% 'square' or 's'  Square
% 'diamond' or 'd' Diamond
% '^' Upward-pointing triangle
% 'v' Downward-pointing triangle
% '>' Right-pointing triangle
% '<' Left-pointing triangle
% 'pentagram' or 'p' Five-pointed star (pentagram)
% 'hexagram' or 'h' Six-pointed star (hexagram)
% 'none' No mark