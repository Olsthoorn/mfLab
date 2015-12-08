function [xUTM,yUTM] = kmlPath2UTM(kmlFName)
%KMLPATH2UTM convert the coordinates in kmlFile to UTM
%
% USAGE: [xUTM,yUTM = kmlPath2UTM(kmlFName)
%
% see also: kmlpath, wgs2utm, utm2wgs
%
% TO 130924

[E   ,N   ] = kmlpath(kmlFName);
[xUTM,yUTM] = wgs2utm(N,E);
