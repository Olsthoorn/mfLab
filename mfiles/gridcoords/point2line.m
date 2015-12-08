function [xpOut,ypOut,dx,dy,dL]=point2line(X,Y,xp,yp)
%POINT2LINE puts point xp yp on the line given by end points of X(1 end) Y(1 end)
%
% Example:
%   [xpOut,ypOut,dx,dy]=point2line(X,Y,xp,yp)
%
%   Useful to shift a point exactly to a line or to compute the distance
%   from the point to the line
%
%   xpOut, ypOut is point on the line
%   dx,dy vector from xp,yp to xpOut,ypOut
%   dL length of vector dx,dy
%
% See also: rotate
%
% TO 100211

L=sqrt((X(end)-X(1)).^2+(Y(end)-Y(1)).^2);
ex=(X(end)-X(1))/L;
ey=(Y(end)-Y(1))/L;


%% intersect drain with line perpendicular to it through Vbak
u=[ex +ey; ey -ex]\[xp-X(1); yp-Y(1)];
xpOut=X(1)+u(1)*ex;
ypOut=Y(1)+u(1)*ey;

dx=xpOut-xp;
dy=ypOut-yp;
dL=sqrt(dx^2+dy^2);
