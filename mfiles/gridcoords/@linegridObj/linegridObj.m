function P = linegridObj(pline,xGr,yGr,zGr,varargin)
%LINEGRIDOBJ --- generates gridLineObj representing the set of line pieces obtained
% by intersecting pline with a grid defined by xGr,yGr,zGr, LAYCBD and minDZ.
%
% USAGE:
%   P = linegridObj(pline,xGr,yGr,zGr[,LAYCBD[,minDZ[,AXIAL]]]);
%
%   where pline is a polyline (1D, 2D or 3D) in the grid specified by the
%   arguments. The number of columns in pline must match the dimension of the
%   grid. so pline=[x] or pline=[x y] or pline=[x y z] (arrays or column vectors
%
%   Returns array of gridLineObj P, in which each element holds the info of the line piece
%   passing through a single cell.
%
%  POut contains
%    x, y, z, ix, iy, iz, idx, xm, ym, zm, L
%
% SEE ALSO: also gridObj/lineObjects gridLineObj cellIndex cellIndices xyzindex inpolygon
%
% TO 100830 130705

[LAYCBD,varargin] = getProp(varargin,'LAYCBD',[]);
if isempty(LAYCBD)
    [LAYCBD,varargin] = getNext(varargin,'double',0);
end

[minDZ,varargin] = getProp(varargin,'minDZ',[]);
if isempty(minDZ)
    [minDZ,varargin] = getNext(varargin,'double',1e-2);
end

[AXIAL,varargin] = getProp(varargin,'AXIAL',[]);
if isempty(AXIAL)
    [AXIAL,~] = getNext(varargin,{'logical','double'},false);
end

gr = gridObj(xGr,yGr,zGr,LAYCBD,minDZ,AXIAL);

P = gr.lineObjects(pline,zRel);
