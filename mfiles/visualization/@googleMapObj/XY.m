function [X Y]=XY(o,E,N)

% [X,Y] = GM.XY(E,N) computes [X,Y] in m relative to center of picture (approximately)
% using the lenght and hight of the lines through the center of the picture
%%
% USAGE:
%   [X Y]=mf_GMpix2m(ix,iy,zoom,px,py)
% 
%   X Y  coordinates in m relative to 0 0 of local tile
%   may be obtained usging
%     [x,y,ix,iy]=mf_GMLL2pix(Lat,Lon,zoomlev);
%%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG
%
% TO 110501

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

if nargin<5, [px py]=ginput; end

[No Eo]=mf_GMpix2LL(ix,iy,zoom,0,0);
fprintf('NoEo=%15.10f %15.10f\n',No,Eo);

n=2^zoom;
w=1/n;

N=(90         -180*w*iy)-py*w*180/256;
E=(360*w*ix-180        )+px*w*360/256;

DN=N-No;  Y=R2*pi*DN/180;
DE=E-Eo;  X=R2*pi*DE/180;

