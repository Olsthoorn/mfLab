function marker=mf_marker(it,markers)
%MF_MARKER yields next in series and switches back to one when length of markers is reached
%
% Example:
%     marker=mf_marker(it[,markers)
%     marker=mf_marker(it,'xo*^>');
%
% if nargin==1, markers='xo*.xsd^v><ph'; white is omitted but can be used by calling
%
% see also: mf_color mf_linetype
%
% TO 120101 120531

markerlist   ='.osphd^v><*x+';

if nargin==1
    NM=length(markerlist);
    ic=rem(it,NM); if ic==0, ic=NM; end
    marker=markerlist(ic);
    return;
else
    if ischar(markers)
        if strcmpi('none',markers)
            markers='w';
        else  % accept colors
            for ic=length(markers):-1:1
                if ~strfind(markerlist,markers(ic)),
                    markers(ic)=[];
                end
            end
            if isempty(markers)
                markers=markerlist;
            end
        end
    elseif isnumeric(markers),
        NM=markers(1);
        if NM<=0,
            markers=markerlist;
        else
            markers=markerlist(1:min(NM,length(markerlist)));
        end
    else
        markers=markerlist;
    end
    NM=length(markers);
    ic=rem(it,NM); if ic==0, ic=NM; end
    marker=markers(ic);
end

    
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