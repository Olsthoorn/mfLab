% TEST_GM2PNG: EXAMPLE of how to use mf_GM2PNG to retrieve figure from Google Maps
%
% USAGE:
%    test_GM2PNG
%
% TO 110501 110514 121008

center='Dubai,AE';      % as name
zoom    = 12;           % between 0 and 21
pixels  = [512 512]';   % each between 1 and 640
maptype = 'roadmap'; % one of

GM = googleMapObj(center,zoom,pixels,maptype,'png');
GM.image;

