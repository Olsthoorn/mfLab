function wYwX=mf_GMpixsizeinm(zoom,N)
%MF_GMPIXSIZEINM computes size of pixel in m according to zoom level
%
% Example:
%   d=mf_GMpixsizeinm(zoom,Northing);
%
%   d is 2 element factor
%   Only Northing is required as extra info
%
%
% SEE ALSO: mf_GMLL2pix mf_GM2PIC mf_GM2PNG mf_GMpix2m
%
% TO 110503

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

n=2^zoom;
w=1/n;

wY=  pi*R2*w/256;
wX=2*pi*R2*w/256*cos(pi*N/180);

wYwX=[wY wX];
