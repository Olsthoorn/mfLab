function L = above(o,varargin)
%BELOW -- Logical L indicting cells with cell mides above plane Z
%
% USAGE:
%    L = gridObj.above('xz',X,Z);
%    L = gridObj.above('yz',Y,Z);
%    L = gridObj,above('xyz',X,Y,Z);
%    L = gridObj,above('z',Z);
%    L = gridObj,above(Z);
%
% code signals what input is provided
%    if xz, then for all y the same line
%    if yz, then for all x the same line
%    if xyz then a true spatial plane is given
%    if z   then Z spans the grid, i.e. has size (gr.Ny,gr.Nx)
%    if code is omitted, same as code == 'z'
%
% SEE ALSO:  gridObj.below, gridObj.intersect
%
% TO 151221

L = ~(o.below(varargin{:}));
