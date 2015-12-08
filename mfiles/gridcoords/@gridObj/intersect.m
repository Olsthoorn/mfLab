function L = intersect(o,varargin)
%BELOW -- find logical L of cells intersected by the given plane Z
%
% USAGE:
%    L = gridObj.intersect('xz',X,Z);
%    L = gridObj.intersect('yz',Y,Z);
%    L = gridObj,intersect('xyz',XY,Z);
%    L = gridObj,intersect('z',Z);
%    L = gridObj,intersect(Z);
%
%
% code signals what input is provided
%    if xz, then for all y the same Z is assumed
%    if yz, then for all x the same Z is assumed
%    if xyz then a true spatial plane is given
%    if z       then Z is a plane spanning the grid, i.e. (gr.Ny,gr.Nx)
%    if omitted then Z is a plane spanning the grid, i.e. (gr.Ny,gr.Nx)
%
% SEE ALSO: gridObj.above, gridObj.below
%
% TO 151122

[code, varargin] = getType(varargin,'char','z');

switch lower(code)
    case 'xz', [X,varargin] = getType(varargin,'double',[]);
               [Z,  ~     ] = getType(varargin,'double',[]);
               X=X(:)'; Z=Z(:)';
               Zm = interp1(X,Z,o.xm,'linear','extrap');
               Zm = bsxfun(@times,Zm,ones(o.Ny,1));
    case 'yz', [Y,varargin] = getType(varargin,'double',[]);  
               [Z,  ~     ] = getType(varargin,'double',[]);
               Y = Y(:); Z=Z(:);
               Zm = interp1(Y,Z,o.ym,'linear','extrap');
               Zm = bsxfun(@times,Zm,ones(1,o.Nx));
    case 'xyz',[X,varargin] = getType(varargin,'double',[]);
               [Y,varargin] = getType(varargin,'double',[]);
               [Z,  ~     ] = getType(varargin,'double',[]);
               Zm = interp2(X(:)',Y(:),Z,o.xm,o.ym,'linear',mean(Z(:)));
    case 'z',  [Zm,  ~     ] = getType(varargin,'double',[]); % skip
    otherwise
        error('code must be one of {''xz'',''yz'',''xyz'',''z'' or omitted} to indicate input');
end

if o.Ny ~= size(Zm,1) || o.Nx ~= size(Zm,2)
   error('size of Z must equal at this point (gr.Ny,gr.Nx) = (%d,%d)',o.Ny,o.Nz);
end

L1 = o.Z(:,:,1:end-1)  > bsxfun(@times,Zm,ones(1,1,o.Nz));
L2 = o.Z(:,:,2:end  ) <= bsxfun(@times,Zm,ones(1,1,o.Nz));
L  = L1 & L2;

