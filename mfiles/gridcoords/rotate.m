function [x,y]=rotate(x,y,x0,y0,alfa)
%ROTATE rotates coordinates over alfa (counter clockwise) degrees around x0,y0
%
% Example:
%   [x,y]=rotate(x,y,x0,y0,alfa)
%
%   TO 090929

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

sx=size(x); sy=size(y); % size x must be equal to size y
  
if sx(1)~=sy(1) || sx(2)~=sy(2)
    error('size of x and y must be same, x0 and y0 scalar, alfa in degrees');
end

alfa=pi*alfa/180;

M=[cos(alfa) sin(alfa); -sin(alfa) cos(alfa)];

xy=[x(:)-x0,y(:)-y0]*M;

x=reshape(xy(:,1)+x0,sx);
y=reshape(xy(:,2)+y0,sy);
