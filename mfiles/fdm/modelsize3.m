function [xGr,yGr,Z,xm,ym,ZM,dx,dy,DZ,Nx,Ny,Nz,LAYCBD]=modelsize3(xGr,yGr,Z,LAYCBD)
%MODELSIZE3 generates grid info grid line coordinates
%
% Rather replace by much more powerfull and versatile gridObj
%
% Example:
%    [xGr,yGr,Z,xm,ym,ZM,dx,dy,DZ,Nx,Ny,Nz,LAYCBD]=modelsize3(xGr,yGr,Z,LAYCBD)
%
% makes values in xGr,yGr,Z unique sorted
% xGr runs upward
% yGr runs downward
% Z   runs downward and is a 3d vector (1,1,Nz)
%
% Output Z is a 3D array with Nx,Ny,Nz+1 elevations of the cell bottoms and
% ceilings.
%
% Input Z may have Nx+1,Ny+1,Nz+1 elevations for all nodes, however the
% output will be as explained before.
%
% If Z input is a vector, Z output will also be a vector. In that case the
% elevations of all cells in a layer will be assumed the same.
%
% if Z input is a vector, so will be DZ and ZM.
%
% Z must have sufficient layer tops and bottoms to match LAYCBD.
%
% Z  output will always be 3D (row, col, layer)
% Z  may be given as
%  1 1D vector in any dimensional direction
%  a 2D array but we wil then assume it is in [nz+1,nx] cross section format
%  a 3D vector or a 3D array with Nx or Nx_1 columns and Ny or Ny+1 rows.
%  The output Z always has size (Ny,Nx,Nz+1).
%  Z in must have the elevation of the cell bottoms and tops of all cells.
%
%  Z output will be ordered from top to bottom and sorted if needed.
%
%  LAYCBD is the (optional) vector NLAY long indicating whether a layer as
%  a confining bed below it. Itis not used here.
%
%
%
% TO 090104 091207 110319 120311

% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%% First make Z 3D if it is not so already

if length(size(Z))<3 && any(size(Z)==1) % then assume Z is given as [nz,1]or [1,nz]
    Z=unique(Z(:));
    Z=Z(end:-1:1);
    Z=reshape(Z,1,1,length(Z));
elseif length(size(Z))==2 % then array is 2D, assume it is a XSection of shape (nz+1,nx)
    Z=permute(Z,[3,2,1]);
else
    % nothing, it is already 3D
end

if xGr(end)    <xGr(1), MUST_FLIP_ALONG_X=1; else MUST_FLIP_ALONG_X=0; end
if yGr(end)    >yGr(1), MUST_FLIP_ALONG_Y=1; else MUST_FLIP_ALONG_Y=0; end
if Z(1,1,end)>Z(1,1,1), MUST_FLIP_ALONG_Z=1; else MUST_FLIP_ALONG_Z=0; end

xGr=unique(xGr(:)');
yGr=unique(yGr(:) ); yGr=yGr(end:-1:1);

if MUST_FLIP_ALONG_X, Z=Z(:,end:-1:1,:); end
if MUST_FLIP_ALONG_Y, Z=Z(end:-1:1,:,:); end
if MUST_FLIP_ALONG_Z, Z=Z(:,:,end:-1:1); end

xm=0.5*(xGr(1:end-1)+xGr(2:end));
ym=0.5*(yGr(1:end-1)+yGr(2:end));

dx=abs(diff(xGr,1,2));
dy=abs(diff(yGr,1,1));

Nx=length(dx);
Ny=length(dy);

if size(yGr,1)~=size(Z,1) && size(yGr,1)~=size(Z,1)+1 && size(Z,1)~=1, error('size(yGr)~=size(Z,1)-1 !'); end
if size(xGr,2)~=size(Z,2) && size(xGr,2)~=size(Z,2)+1 && size(Z,2)~=1, error('size(xGr)~=size(Z,2)-1 !'); end

% Now look in xGr and yGr direction
if size(Z,1)==size(yGr,1)
    Z=0.5*(Z(1:end-1,:,:)+Z(2:end,:,:));
end
if size(Z,2)==size(xGr,2)
    Z=0.5*(Z(:,1:end-1,:)+Z(:,2:end,:));
end

ZM=0.5*(Z(:,:,1:end-1)+Z(:,:,2:end));
DZ=abs(diff(Z,1,3));

%% Verify LAYCBD compatibility
if nargin<4,
    LAYCBD=zeros(size(Z(1:end-1)));
end
[isAqf,LAYCBD]=isAquifer(size(Z,3)-1,LAYCBD);

Nz=sum(isAqf>0);

ZM=ZM(:,:,isAqf>0);


