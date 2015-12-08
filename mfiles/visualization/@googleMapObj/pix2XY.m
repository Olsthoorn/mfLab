function [X Y]=mf_GMpix2XY(ix,iy,zoom,px,py)

% GMpix2XY: computes dist in m to zero meridian and equator
%
% USAGE:
%   [X Y]=mf_GMpix2LL(ix,iy,zoomlev,LATC,LONC,xpix,ypix);
%   [X Y]=mf_GMpix2LL(ix,iy,zoomlev,LATC,LONC);    % invokes ginput to get xpix ypix
% 
%   X,Y  distance along face of earth to zero meridian and equator
%   ix=GM x-tile number
%   iy=GM y-yile number
%%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110501

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

if nargin<5, [px,py]=ginput; end

n=2^zoom;
w=1/n;

y=w*(iy+py/256);
x=w*(ix+px/256);

N=mf_GMy2lat(y);
E=mf_GMx2lon(x);

Y=pi/180*N*R2;
X=pi/180*E*R2*cos(pi/180*N);

fprintf('%12.2f %12.2f\n',y,x);
fprintf('%12.2f %12.2f\n',N,E);
fprintf('%12.2f %12.2f\n',Y,X);
fprintf('\n');

