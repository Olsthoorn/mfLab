function [A,xGr,yGr,Ix,Iy,meta]=NHI_readASC(fname,xLim,yLim)
%NHI_READASC reads NHI ASCII datafile (ESRI), select between given coordinates
%
% Example:
%    [A,xGr,yGr,xm,ym,Ix,Iy]=readASC(FName[,pth[,xLim,yLim [,meta]]])
%    [A     ,xGr,yGr,xm,ym,Ix,Iy]=NHI_readASC(FName,pth,xLim,yLim)
%    [meta,xGr,yGr,xm,ym,Ix,Iy]=NHI_readASC(FName,pth,xLim,yLim,'meta')
%
%FILE STRUCTURE
%   Structure
%   ncols         1200
%   nrows         1300
%   xllcorner     0                sometimes xllcenter
%   yllcorner     300000           sometimes yllcenter
%   cellsize      250
%   NODATA_value  -9999
%
% See also: NHI_read getNHImeta getNHIBCN getNHIASC
%
%  TO 2009 110425


fprintf('Reading file    ''%s''\n',fname);

fid=fopen(fname,'r');  if fid<0, error('can''t open file <<%s>>',fname); end

header = textscan(fid, '%s %f', 6);

for i=1:size(header{1},1)
    meta.(header{1}{i})=header{2}(i);
end

%% Set xGr and yGr of cell centers
xGr=meta.XLLCORNER+meta.CELLSIZE*(0:(meta.NCOLS));                    
yGr=meta.YLLCORNER+meta.CELLSIZE*(0:(meta.NROWS));

xm=0.5*(xGr(1:end-1)+xGr(2:end));
ym=0.5*(yGr(1:end-1)+yGr(2:end));

A=fscanf(fid,'%f',[meta.NCOLS,meta.NROWS])';
A(A==meta.NODATA_VALUE)=NaN;

if nargin<3
    xLim=meta.XLLCORNER+[0, meta.NCOLS*meta.CELLSIZE];
    yLim=meta.YLLCORNER+[0, meta.NROWS*meta.CELLSIZE];
end

Ix=find(xm>xLim(1) & xm<xLim(end));
Iy=find(ym>yLim(1) & ym<yLim(end));

A=A(Iy,Ix);

fclose(fid);

