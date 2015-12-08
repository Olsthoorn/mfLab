function [subArea,e,n]=subArea(Area,LONLIM,LATLIM,lonlim,latlim)
%SUBAREA cuts subarea from area given LONLIM LATLIM of area and lonlim latlim of subArea
%
% Example:
%    [subArea,e,n]=subArea(Area,LONLIM,LATLIM,lonlim,latlim)
%
% See also: mf_getdemfromtiff
%
% TO 110523

nx=size(Area,2); dx=diff(LONLIM)/nx;
ny=size(Area,1); dy=diff(LATLIM)/ny;

Ix=round((lonlim-LONLIM(1))/diff(LONLIM)*nx); Ix=Ix(1):Ix(2);
Iy=round((latlim-LATLIM(1))/diff(LATLIM)*ny); Iy=Iy(1):Iy(2);

e=LONLIM(1)+(Ix(:)'-1)*dx;
n=LATLIM(1)+(Iy(:) -1)*dy;


subArea=Area(Iy,Ix);

