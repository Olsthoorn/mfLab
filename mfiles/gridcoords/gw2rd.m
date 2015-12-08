function [RDx,RDy]=gw2rd(GWx,GWy)
%GW2RD convert local Amsterdam Water Supply Dune Area coords to Dutch national coordinates.
%
% gw2rd.m: transforms gw-coordinates into rd-coordinates
%               
% Example
%    [RDx RDy]=gw2rd(GWx,GWy)
%
% See also: rd2gw wgs2rd rd2wgs
%
% Copyright P.Kamps 10/nov/1999


cosa = 0.923319;
sina = 0.384034;
%xgw  = 61892.94;
%ygw  = 885.364;
xca  = 155000;
yca  = 463000;
xrn  = -57486.9;
yrn  = 22951.5;

RDx  = xrn + (GWx * cosa ) + (GWy * sina ) + xca;
RDy  = yrn + (GWy * cosa ) - (GWx * sina ) + yca;

