function meta=getNHImeta(fid)
%GETNHIMETA gets the 6 lines of meta data from any of the NHI .ASC files.
%
% Example:
%    meta = get(NHImeta(fid)  -- meta data from file ...fname.ASC ...
%
% Meta data (first 6 line sof ASCI file)
%   Structure
%   ncols         1200
%   nrows         1300
%   xllcorner     0                sometimes xllcenter
%   yllcorner     300000           sometimes yllcenter
%   cellsize      250
%   NODATA_value  -9999
%
% EXAMPLE:
%   meta=NHI_readASC(fname)
%   [meta,xGr,yGr]=NHI_readASC(fname)
%   [meta,xGr,yGr,Ix,Iy]=NHI_readASC(fname,xLim,yLim)
%
%  TO 2009 110425

header = textscan(fid, '%s %f', 6);

for i=1:size(header{1},1)
    meta.(header{1}{i})=header{2}(i);
end

