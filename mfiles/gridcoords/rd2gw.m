function [GWx,GWy]=rd2gw(RDx,RDy)
%RD2GW converts Duth national coords to local coords of Amsterdam Water Supply Dune Area
%
% Example:
%    [GWx GWy]=RD2GW(RDx,RDy)
%
% See also: gw2rd wgs2rd rd2wgs
%
% Copyright P.Kamps 27/okt/1999

cosa = 0.923319;
sina = 0.384034;
xgw  = 61892.94;
ygw  = 885.364;
xca  = 155000;
yca  = 463000;

GWx  = xgw + ( RDx - xca ) * cosa - ( RDy - yca ) * sina;
GWy  = ygw + ( RDy - yca ) * cosa + ( RDx - xca ) * sina;
