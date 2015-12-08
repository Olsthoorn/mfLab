function hdl=mf_logo(varargin)
%MF_LOGO plots string www.google.com/p/mflab on lower left of figure
%
% USAGE:
%    hdl=mf_logo(ax,xOff,yOff,varargin)
%
% Example:
%   [hdl]=mf_logo;                % generate logo give handle to it
%   [hdl]=mf_logo(ax);            % same on axis ax
%   set(hdl,'string',text);       % put this text on the logo
%   set(hdl,'position',[xOff,yOff]);     % put logo on this location
%   [hdl]=mf_logo(...,property,value,property,value,...);
%
% Example:
%   [hdl]=mf_logo( ax,xOff,yOff,'units','normalized');
%   [hdl]=mf_logo( ax,xOff,yOff,'units','data');
%
%  See also
%   Properties of graphic object 'text' for list of properties.
% 
% TO 110422  120409

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

s='www.code.google.com/p/mfLab';

[ax       ,varargin] = getNext(varargin,'axis' ,  gca);
[xOff     ,varargin] = getNext(varargin,'double',0.05);
[yOff     ,varargin] = getNext(varargin,'double',xOff);

[units    ,varargin] = getProp(varargin,'units','normalized');
[fontname ,varargin] = getProp(varargin,'fontname','courier');
[fontangle,varargin] = getProp(varargin,'fontAngle','italic');
[color    ,varargin] = getProp(varargin,'color','k');

defaultProps={'units',units,'fontname',fontname,'fontangle',fontangle,'color',color};

hdl=text(xOff,yOff,s,defaultProps{:},varargin{:},'parent',ax);

end

%% Legal properties of text
% 	BackgroundColor
% 	Color
% 	EdgeColor
%	Extent = [-1.93237 8.55615 21.9002 2.49554]
% 	FontAngle: [ {normal} | italic | oblique ]
% 	FontName
% 	FontSize
% 	FontUnits: [ inches | centimeters | normalized | {points} | pixels ]
% 	FontWeight: [ light | {normal} | demi | bold ]
% 	HorizontalAlignment: [ {left} | center | right ]
% 	LineStyle: [ {-} | -- | : | -. | none ]
% 	LineWidth
% 	Margin
% 	Position
% 	Rotation
% 	String
% 	Units: [ inches | centimeters | normalized | points | pixels | characters | {data} ]
% 	Interpreter: [ latex | {tex} | none ]
% 	VerticalAlignment: [ top | cap | {middle} | baseline | bottom ]
% 	Parent
% 	Tag
%   Type=text
% 	UserData
% 	Visible: [ {on} | off ]
