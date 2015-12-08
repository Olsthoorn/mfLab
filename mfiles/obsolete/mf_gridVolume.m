function DV=mf_gridVolume(xGr,yGr,zGr,AXIAL)
%MF_GRIDVOLUME compute the cell volumes of the finite difference grid
%  assuming xGr,yGr and zGr are vectors.
%  AXIAL has value 0 or non zero
% USAGE
%   DV= mf_gridVolume(xGr,yGr,zGr[,AXIAL]);
%
% TO 120104
%

if nargin<4, AXIAL=0; end

[~,~,~,xm,~,~,Dx,Dy,Dz,NX,NY,NZ]=modelsize3(xGr,yGr,zGr);

if AXIAL
    DY=repmat(2*pi*xm,[NY,1,NZ]);
else
    DY=repmat(Dy,[1,NX,NZ]);
end
DX=repmat( Dx    ,[NY,1,NZ]);
DZ=repmat( Dz    ,[NY,NX,1]);
DV    =DX.*DY.*DZ;
