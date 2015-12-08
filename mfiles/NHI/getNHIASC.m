function [A,meta]=getNHIASC(fname,Ix,Iy)
%GETNHIASC reads ASCII (ESRI) datafile, select between given coordinates
%
% Example:
%   [A,meta]=getNHIASC(fname,Ix,Iy)
%
% TO 120427

fprintf('Reading file    ''%s''\n',fname);

fid=fopen(fname,'r');  if fid<0, error('can''t open file <<%s>>',fname); end

header = textscan(fid, '%s %f', 6);

 for i=1:size(header{1},1)
     meta.(header{1}{i})=header{2}(i);
 end

%% Set xGr and yGr of cell centers

% The first value in the file must by the top-left column with indices Ix=1
% and Iy=1 as in the MODFLOW grid. This is the case with the ASCII input
% files of the NHI. So no flipping is necessary
A=fscanf(fid,'%f',[meta.NCOLS,meta.NROWS])';
A(A==meta.NODATA_VALUE)=NaN;

A=A(Iy,Ix);

fclose(fid);
