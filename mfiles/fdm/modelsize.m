function [x,y,xm,ym,dx,dy,Nx,Ny]=modelsize(x,y)
%MODELSIZE generates grid info from gridLine coordinates
%
%[x,y,xm,ym,dx,dy,Nx.Ny]=modelsize(x,y)
%
% What does it do?
%    makes values in x and y unique and sorted with correct orientation
%    then compute the other values
%
% xm,ym cell center cordinates
% dx, dy cell widths (positivie)
% Nx, Ny size grid 
%
% See also: gridObj modelsize3 fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% replace by gridObj
%
% TO 2004 120804

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

x=unique(x(:)');
if isvector(y)
    y=sort(unique(y(:)),'descend');
else
    if y(1,end)>y(1,1)
        y = flipud(y);
    end
end
xm=0.5*(x(1:end-1)  +x(2:end)  );
ym=0.5*(y(1:end-1,:)+y(2:end,:));

dx=    diff(x,1,2);
dy=abs(diff(y,1,1));

Nx=size(dx,2);
Ny=size(dy,1);
