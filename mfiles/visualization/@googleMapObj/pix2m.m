function [Lat Lon]=pix2m(o,px,py)
% GMpix2LL: computes Lat Lon from local pix coordinates
%%
% USAGE:
%   [Lat Lon]= GM.pix2m(px,py)
% 
% px,py are local pixell coordinates in tile
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG
%
% TO 110501 12108

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

if nargin<5, [px py]=ginput; end

[No Eo]=o.pix2LL(0,0);
fprintf('NoEo=%15.10f %15.10f\n',No,Eo);

n=2^zoom;
w=1/n;

N=(90         -180*w*iy)-py*w*180/256;
E=(360*w*ix-180        )+px*w*360/256;

DN=N-No;  Y=R2*pi*DN/180;
DE=E-Eo;  X=R2*pi*DE/180;

