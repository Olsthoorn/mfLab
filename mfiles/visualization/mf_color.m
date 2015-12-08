function clr=mf_color(it,clrs)
%MF_COLOR yields next in color series and restarts after last color has been used
%
% Example:
%    clr=mf_clr(it);            % uses default colors
%    clr=mf_clr(it,clrs)        % uses your colors
%    clr=mf_clr(it,'brgkwyb');  % uses specified colors
%    default colors: clrs='brgkmcy'
%    white is omitted but can be used by calling
%
% See also: mf_linetype mf_marker
%
% TO 120101 120531

clrlist      = 'brgkmcy';

if nargin==0, clr=clrlist; return; end

if nargin==1
    NC=length(clrlist);
    ic=rem(it,NC); ic(ic==0)=NC;
    clr=clrlist(ic);
else 
    if ischar(clrs)
        if strcmpi('none',clrs)
            clrs='w';
        else  % accept colors
            for ic=length(clrs):-1:1
                if ~strfind(clrlist,clrs(ic)),
                    clrs(ic)=[];
                end
            end
            if isempty(clrs)
                clrs=clrlist;
            end
        end
    elseif isnumeric(clrs),
        NC=clrs(1);
        if NC<=0,
            clrs=clrlist;
        else
            clrs=clrlist(1:min(NC,length(clrlist)));
        end
    else
        clrs=clrlist;
    end
    NC=length(clrs);
    ic=rem(it,NC); if ic==0, ic=NC; end
    clr=clrs(ic);
end
            
% Matlab's linespecs:    
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