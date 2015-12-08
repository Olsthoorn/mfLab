% TEST_GM2PNG: EXAMPLE of how to use mf_GM2PNG to retrieve figure from Google Maps
%
% USAGE:
%    test_GM2PNG
%
% TO 110501 110514 121008

%% We start making a struct holding the parameter values that Goolge Maps

center='Brooklyn+Bridge,New+York,NY';
zoom    = 13;           % between 0 and 21
pixels  = [640 640]';   % each between 1 and 640
maptype = 'hybrid'; % one of

GM = googleMapObj(center,zoom,pixels,maptype,'png');
GM.image;

