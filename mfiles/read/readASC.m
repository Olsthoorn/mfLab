function [A,xGr,yGr,xm,ym,Ix,Iy]=readASC(FName,xLim,yLim)
%%READASC -- reads ESRI ascii files (georeferenced by meta data)
%
% USAGE:
%   [A,xGr,yGr,xm,ym,IxLim,IyLim]=readASC(FName,xLim,yLim [,header])
%
% Stucture of ESRI ASCII file: first 6 lines contain keyword and value
% ncols         value
% nrows         value
% xllcorner | xllcenter    value
% yllcorner | yllcenter    value
% cellsize      value
% NODATA_value  value
%
% TO 010109

%% Open
fprintf('Reading file    ''%s''\n',FName);

try
    fid=fopen(FName,'r');
    if fid<1
        error('%s: Can''t open file %s',FName);
    end
catch Me
    fid = fopen([FName,'asc'],'r');
    if fid<1
        error([Me.message,'Can''t open file <<%s>>'],FName);
    end
end
    
%% Assert there are no decimal comma's
if any(fscanf(fid,'%c',4000)==',')
    error('%s: file %s contains decimal comma''s instead of periods',mfilename,FName);
end

%% Reset file pointer
fseek(fid,0,-1);

%% Read meta data
keyw = cell(1,6);
val  = NaN(1,6);

for i=1:6
    keyw{i} = fscanf(fid,'%s',1);
    val(i)  = fscanf(fid,'%f',1);
end

ncols    = val(strmatchi('ncols'       ,keyw));
nrows    = val(strmatchi('nrows'       ,keyw));
xll      = val(strmatchi('xll'         ,keyw));
yll      = val(strmatchi('yll'         ,keyw));
cellsize = val(strmatchi('cellsize'    ,keyw));
NODATA   = val(strmatchi('NODATA_value',keyw));

%% Generate grid coordinates

if strmatchi('xllcorner',keyw)
    xGr=[xll, xll+cellsize*(1:ncols)];
else % xllcenter
    xGr=[xll, xll+cellsize*(1:ncols)] -0.5*cellsize;
end

if strmatchi('yllcorner',keyw)
    yGr=[yll, yll+cellsize*(1:nrows)];
else % yllcenter
    yGr=[yll, yll+cellsize*(1:nrows)] - 0.5*cellsize;
end

yGr=yGr(end:-1:1);  % bacause top lines has highest y coordinate

xm=0.5*(xGr(1:end-1)+xGr(2:end));
ym=0.5*(yGr(1:end-1)+yGr(2:end));

%% Read and process file
A=fscanf(fid,'%f',[ncols,nrows])';
A(A==NODATA)=NaN;

%% Default xLim and yLim if not specified in call
if nargin<3, xLim=xGr([1 end]); yLim=sort(yGr([1 end])); end

%% Extract desired portion from file
Ix=find(xm>xLim(1) & xm<xLim(2));
Iy=find(ym>yLim(1) & ym<yLim(2));

xGr=xGr(xGr>=xLim(1) & xGr<=xLim(2));
yGr=yGr(yGr>=yLim(1) & yGr<=yLim(2));

xm=xm(Ix);
ym=ym(Iy);

A=A(Iy,Ix);

fclose(fid);

