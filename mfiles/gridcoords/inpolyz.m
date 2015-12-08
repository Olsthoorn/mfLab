function [IN,ON] = inpolyz(X,Z,xv,zv)
%INPOLYZ returns logical arrray for zx plane telling which cells are in vertical polygon
%
% Example:
%    IN = inpolyz(X,Y,XV,YV) returns a matrix IN the size of X and Z..
%    IN(p,q) = 1 if the point (X(p,q), Z(p,q)) is either strictly inside
%
% See also: inpolygon inpoly
%
% TO 130409

X = XS(X);
Z = XS(Z);

[IN,ON] = inpolygon(X,Z,xv,zv);

IN = XS(IN);
ON = XS(ON);