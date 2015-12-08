function o = head(o,varargin)
% function animateObj.head(gr,well,,varargin);
% plot head XZ YZ or XY plane, use property value pairs to specify which plane
% varargin is propty value pair showing what to simulate
% well will be plotted unless omitted
%
% exmaples
%  'ix',ix  ---> contour column iz
%  'iy',it  ---> contour row    iy
%, 'iz',iz  ---> contour plane  iz
%
% TO 121219

error(['%s: this function is obsolete\n',...
    'REMEDY:  use one of these;\n',...
    '         animateObj.headXS  (plots an X-section, i.e. along the x-axis)\n',...
    '         animateObj.headYS  (plots an Y-section, i.e. along the y-axis)\n',...
    '         animateObj.headXY  (plots a  Z-section, or a XY slices, i.e. horizontal)\n'...
    '         TO 130403'],  mfilename);
    
