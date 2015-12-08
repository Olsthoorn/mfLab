function h=colorCell(I,x,y,clr)
%COLORCELL fills cells with color
%
% USAGE:
%    h=colorCell(I,x,y,clr)
%
% cells specified by number global cell nr I and grid lines x and y will be colored
% using the the given colorSpec clr.
%
% TODO: seems outdated, nowadays I would use the gridObj TO 130427
%
% TO 090315

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


dx=diff(x);
dy=diff(y);
Ny=length(dy);

for i=1:length(I)
    iy=rem(I(i),Ny); if iy==0; iy=Ny; end
    ix=(I(i)-iy)/Ny+1;
    if dy(iy)>0
        h=rectangle('position',[x(ix),y(iy  ),dx(ix), dy(iy)],...
            'edgecolor',clr,'facecolor',clr);
    else
        h=rectangle('position',[x(ix),y(iy+1),dx(ix),-dy(iy)],...
            'edgecolor',clr,'facecolor',clr);
    end
end
