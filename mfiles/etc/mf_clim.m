function [c1c4]=mf_clim(c2,c3,i2,i3,L)
%MF_CLIM yields the extremes of the colormap values to be used as clim([c1 c4])
%
% USAGE:
%    [c1 c4]=mf_clim(c2,c3,i2,i3,L)
% 
% Serves to plot the object using only the colors in i2...i3 of the
% colormap which has values c2,c3.
%
% will be obsolete in the future
% 
% See also setmulticolormap
%
% TO 110320

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

i1=1; i4=L;
c1c4=[i4-i2 i2-i1; i4-i3 i3-i1]/(i4-i1)\[c2; c3];
