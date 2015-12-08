function Var=rarray(fid,RowCol,norec)
%RARRAY general MODFLOW/MT3DMS array reader
%
% Example:
%     Var=rarray(fid,RowCol,norec)
%
% we will use this rarray in the future as the general array reader for
% MODFLOW, MT3DMS, SEAWAT and compatible models
% I'm not sure if it works right now. It's writing equivalent is rarray,
% which is the general array writer to produce input files for MODFLOW etc.
%
% See also: warray

% TO 110112

if ~exist('norec','var') % || isempty(norec)
    Var=mudread(fid,RowCol);
else
    Var=mudread(fid,RowCol,norec);
end
