function [N E]=mf_GMpix2LL(ix,iy,zoom,px,py)
%MF_GMPIX2LL computes lat long from local pix coordinates given GM tile and zoom level.
%
% Example:
%   [LAT LON]=mf_GMpix2LL(ix,iy,zoomlev,LATC,LONC,xpix,ypix);
%   [LAT LON]=mf_GMpix2LL(ix,iy,zoomlev,LATC,LONC);    % invokes ginput to get xpix ypix
% 
%   LAT,LONC  are lat and lon of center of picture
%   ix=GM x-tile number
%   iy=GM y-yile number
%%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110501

if nargin<5, [px,py]=ginput; end

n=2^zoom;
w=1/n;

y=w*(iy+py/256);
x=w*(ix+px/256);

N=mf_GMy2lat(y);
E=mf_GMx2lon(x);
