function showshp(shpfile)
%SHOWSHP plots shapfile in black by reading appropriate shapefiles
%
% Example:
%    showshp(shpfile)
%
% needs readshp and dbfread must be in the path
% the <shapfile>.shp, <shapfile>.shx, <shapfile>.dbf  must be available
%
% See alo: readshp
%
% TO 090722

shape=readshp([shppath shpfile]);

for i=1:length(shape)
    plot(shape(i).x{1},shape(i).y{1},'k');
    hold on
end