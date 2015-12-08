function [px,py,ix,iy]=mf_GMLL2pix(N,E,zoom)
%GMLL2PIX translates lat lon to pix coordinates on local GM tile
%
% Example:
%   [px,py,ix,iy]=mf_GMLL2pix(N,E,zoom)
%
%   input:N= LAT E=LON and the used zoom leve 0-21
%   px,py are pixel coordinates from the LAT LON
%   to the center of the tile that contains this
%   point at the given zoom level.
%
%   ix,iy is the tile level/position in Goolgle Maps tile system
%   zoom = 0-21 according to google
%
% SEE ALSO: mf_GM2PIC mf_GMpix2LL
%
% TO 110501

n=2^zoom;         %fprintf('%15.10g',n);
w=1/n;            %fprintf('%15.10g',w)

y=mf_GMlat2y(N);
x=mf_GMlon2x(E);

ix=fix(x/w);
iy=fix(y/w);

px =256*(x/w-ix);  %fprintf('%15.10g',px);
py =256*(y/w-iy);  %fprintf('%15.10g',py);
                  %fprintf('\n'); 