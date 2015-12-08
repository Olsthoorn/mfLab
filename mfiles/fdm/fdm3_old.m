function [Phi,Q,Qx,Qy,Qz]=fdm3(x,y,Z,kx,ky,kz,IBOUND,IH,FQ) %% future: ,hdr,Cdr,hrv,hbr,Crv,hgh,Cgh)
%FDM3 a 3D block-centred steady-state finite difference model
%
% Example:
%    [Phi,Q,Qx,Qy,Qz]=fdm3(x,y,Z,x,kx,ky,kz,IBOUND,IH,FQ)
%
% INPUT:
%    x,y,Z mesh coordinates
%    Z may be a vector, a 3D vertical vector or a full size 3D array of size Ny,Nx,Nz+1
%    kx,ky,kz conductivity arrays
%    IBOUND as in MODFLOW
%    IH=initial heads (NaN for ordinary points),
%    FQ=fixed nodal flows
% OUTPUT:
%    Phi,Q computed heads and cell balances
%    Qx,Qy,Qz is cell face flow, positive along positive axis direction
%    Future:  hdr=elevation drains, Cdr=resistance drains
%             hrv=elevation river , Crv=resistance river, hbr=elevation bottom river
%             hgh=elevation general head, Cgh is resistance general head
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 070513
% TO 080226 implemented inactive cells and true fixed heads
% TO 090216 small edits

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

gr = gridObj(x,y,Z);

Nodes = reshape(1:gr.Nxyz,gr.size);               % Node numbering
IW=Nodes(:,1:end-1,:);    IE=Nodes(:,2:end,:);
IN=Nodes(1:end-1,:,:);    IS=Nodes(2:end,:,:);
IT=Nodes(:,:,1:end-1);    IB=Nodes(:,:,2:end);

% resistances and conducctances
warning('off','all');  % to allow zero conductivities for inactive cells
RX=0.5*gr.DX./(gr.DY.*gr.DZ.*kx); Cx=1./(RX(:,1:end-1,:)+RX(:,2:end,:));
RY=0.5*gr.DY./(gr.DZ.*gr.DX.*ky); Cy=1./(RY(1:end-1,:,:)+RY(2:end,:,:));
RZ=0.5*gr.DZ./(gr.DX.*gr.DY.*kz); Cz=1./(RZ(:,:,1:end-1)+RZ(:,:,2:end));
warning('on','all');

A=sparse([IE(:);IW(:);IN(:);IS(:);IT(:);IB(:)],...
         [IW(:);IE(:);IS(:);IN(:);IB(:);IT(:)],...
        -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],...
         gr.Nxyz,gr.Nxyz,7*gr.Nxyz);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct = IBOUND~=0; % active cells
I    = IBOUND >0; % active cells but not fixed heads = cells with heads to be computed
Ifh  = IBOUND <0; % active cells with fixed heads


Phi=IH(:); FQ =FQ(:); Q  =FQ;

if any(I(:))  % something to compute
    Phi( I)=spdiags(Adiag(I   ),0,A(I   ,I   ))\(FQ(I)-A(I,Ifh)*Phi(Ifh)); % solve
end
if any(IAct)
    Q(IAct)=spdiags(Adiag(IAct),0,A(IAct,IAct))*Phi(IAct );		           % reshape
end

Phi= reshape( Phi,size(IBOUND));
Q  = reshape( Q  ,size(IBOUND));

if nargout>2
    Qx=-Cx.*diff(Phi,1,2)*sign(gr.xGr(end)-gr.xGr(1)); Qx(isnan(Qx))=0;  % Flow across horizontal cell faces
    Qy=-Cy.*diff(Phi,1,1)*sign(gr.yGr(end)-gr.yGr(1)); Qy(isnan(Qy))=0;  % Flow across vertical cell faces
    Qz=-Cz.*diff(Phi,1,3)*sign(gr.Z(end)-gr.Z(1));     Qz(isnan(Qz))=0;  % Flow across vertical cell faces
end