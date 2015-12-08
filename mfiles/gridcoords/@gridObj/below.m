function L = below(o,varargin)
%BELOW -- find logical array L indicating cells below a plane given by [X[,Y],Z
%
% USAGE:
%    L = gridObj.below('xz',x,z); % x,z are vectors
%    L = gridObj.below('yz',y,z); % y,z are vectors
%    L = gridObj.below('xyz',x,y,z); % x,y vectors, z array spanning x, and y
%    L = gridObj.below('z',z); % z is a plane array spanning x, and y
%    L = gridObj.below(z);     % z is array spanning gridObj.xm, gridObj.ym 
%
% code signals what input is provided
%    if xz, then for all y the same line
%    if yz, then for all x the same line
%    if xyz then a true spatial plane is given
%    if z        then z spans the grid implied by the gridObj
%    if omitted, then z spans the grid implied by the gridObj.
%
% SEE ALSO: gridObj.above, gridObj.intersect
%
% TO 151121

[code,varargin] = getType(varargin,'char','z');

switch lower(code)
    case 'xz', [X,varargin] = getType(varargin,'double',[]);
               [Z,  ~     ] = getType(varargin,'double',[]);
               X=X(:)'; Z=Z(:)';
               Zm = interp1(X,Z,o.xm,'lineair','extrap');
               Zm = bsxfun(@times,Zm,ones(o.Ny,1));
    case 'yz', [Y,varargin] = getType(varargin,'double',[]);  
               [Z,  ~     ] = getType(varargin,'double',[]);
               Y = Y(:); Z=Z(:);
               Zm = interp1(Y,Z,o.ym,'lineari','extrap');
               Zm = bsxfun(@times,Zm,ones(1,o.Nx));
    case 'xyz',[X,varargin] = getType(varargin,'double',[]);
               [Y,varargin] = getType(varargin,'double',[]);
               [Z,  ~     ] = getType(varargin,'double',[]);
               Zm = interp2(X(:),Y(:)',Z,o.xm,o.ym,'linear',mean(Z(:)));
    case 'z',  [Zm,  ~     ] = getType(varargin,'double',[]);
                if ~all(size(Zm) == o.size)
                    error('provided Z must be of size of grid, i.e. (gr.Ny,gr.Nx) = (%d,%d)',...
                        o.Ny,o.nx);
                end
    otherwise
        error('code must be one of {''xz'',''yz'',''xyz'' or ''z''} to indicate input');
end

L = o.ZM < bsxfun(@times,Zm,ones(1,1,o.Nz));

