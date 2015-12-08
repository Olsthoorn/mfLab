function hdl=gridsurf(xGr,yGr,ZM,C)
%GRIDSURF plots a surface using grid coordinates for x and y and center of cell values for ZM
%
% Example:
%    hdl=gridsurf(xGr,yGr,ZM)
%
%   A bit special, normally the coordates and values are all given at cell
%   centers. In this case we first iterpolate the center values to the
%   corners and the proceed in contouring the values at the corners.
%
%   An alternative is using gridObj.xp and gridObj.xp, gridObj.zp
%
%  See also: gridObj
%
% TO 091225

[ny,nx]=size(ZM);

if nx~=size(xGr,2)-1 || ny~=size(yGr,1)-1
    error('size of ZM, (nx=%d, ny=%d) does not match length of xGr-1=%d or length of yGr-1=%d!\n',...
        nx,ny,size(xGr,2)-1,size(yGr,1)-1);
end

switch nargin
    case 4, hdl=surf(xGr,yGr,mid2corner(ZM),C);
    case 3, hdl=surf(xGr,yGr,mid2corner(ZM),ZM);
    otherwise
        error('gridsurf: needs at least 3 inputs not %d',nargin);
end

function Z=mid2corner(ZM)
% Z=mid2corner(ZM)
% Computes the Z of corner points of a grid of which the Z values of the
% centers are given. It does so approximately by taking the average of the
% surrounding ZM around each grid point in Z
% necessary to plot a surface with given xGr, yGr ZM where ZM has
% one column and row less than xGr and yGr as is always the case in a
% finite element model.

[ny,nx]=size(ZM);

 Z=zeros(ny+1,nx+1);
 C =zeros(ny+1,nx+1);
 for i=[0 1]
     for j=[0 1]
         Z((1:ny)+i,(1:nx)+j)=Z((1:ny)+i,(1:nx)+j)+ZM;
         C((1:ny)+i,(1:nx)+j)=C((1:ny)+i,(1:nx)+j)+ 1;
     end
 end
 Z=Z./C;
 

