function bound=digitiz(z,color)
%DIGITIZ digitize on screen axis using mouse and show what has been selected
%
% USAGE:
%    bound=digitiz(z,color)
%
%    Using another button finishes this function and yields
%    a 3 column matrix with [x y z]
%
%    During digitizing, the drawn line is visible in the given color
%
% TO 070101

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


button=1; x=[]; y=[];
while button==1
    [xp,yp,button]=ginput(1);
    x=[x;xp]; y=[y;yp];
    if length(x)==1,
        h=line(x,y);
    else
        set(h,'xdata',x,'ydata',y,'color',color); drawnow;
    end
end
bound=[x,y,ones(length(x),1)*z];
