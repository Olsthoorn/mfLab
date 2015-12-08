function [xRD,yRD]=kmlpath2rd(kmlfname)
%KMLPATH2RD give xRD and yRD coordinates of GE path in kml file
%
% Example:
%   [xRD,yRD]=kmlpath2rd(kmlfname)
%
% where
%   kmlfname = fname given by you to your savd as kmfl file Google Earth path
%   [xRD,yRD] "Rijksdriehoek coordinates" Dutch national x,y system
%
% Example (recipie):
%   1) launch GE, select add path, click a path
%   2) when finished save it as a kml file in a convenient directory
%      e.g. as "mypath.kml"
%   3) type
%      [xRD,yRD]=kmlpath2rd('mypath');
%
% See also kmlpath wgs3rd rd2wgs
%
% TO 091215 110427

[E,N]=kmlpath(kmlfname);

[xRD,yRD]=wgs2rd(E,N);

