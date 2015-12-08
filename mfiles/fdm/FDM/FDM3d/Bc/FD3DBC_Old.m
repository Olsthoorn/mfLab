function [Phi,Q,QRF,QFF,QLF]=FD3DBC(x,y,z,kx,ky,kz,FH,Q)
% function [Phi,Q,QRF,QFF,QLF]=D3DBC(x,y,z,kx,ky,kz,FH,Q)
% Defines and solves a 3D block-centred finite difference model (MODFLOW-type model)
% x,y,z mesh coordinates (model cell side coordinates)
% kx,ky,kz cell conductivities,
% FH=fixed heads for cells (NaN for ordinary points)
% Q=fixed cell flows, injectionsn positive, extractions negative, else zero
% Phi is computed heads for cells
% Q computed cell balances
% QRF=Qright face, QFF=Q front face, QLF=Q lower face zoals in MODFLOW
% TO 991119, TO 000120 TO 010517

HUGE=1e20;

x=x(:)';	Nx=length(x)-1;	dx=abs(diff(x));
y=y(:);	Ny=length(y)-1; 	dy=abs(diff(y));
z=z(:)';	Nz=length(z)-1;	dz=abs(diff(z)); Ncel=Nx*Ny*Nz;
[DX,DY,DZ]=meshgrid(dx,dy,dz);

if isempty(FH),  FH=NaN*zeros(Ny,Nx,Nz); end; FH=FH(:);
if isempty( Q),   Q=    zeros(Ny,Nx,Nz); end;  Q= Q(:);

% node numbering
Nodes = reshape([1:Nx*Ny*Nz],Ny,Nx,Nz);
Ix=Nodes(:,1:end-1,:); Jx=Nodes(:,2:end  ,:);
Iy=Nodes(1:end-1,:,:); Jy=Nodes(2:end  ,:,:);
Iz=Nodes(:,:,1:end-1); Jz=Nodes(:,:,2:end  );

Rx=0.5*DX./(DY.*DZ)./kx;			% resistance of flow for half cell length [d/m2]
Ry=0.5*DY./(DZ.*DX)./ky;			% same for y direction
Rz=0.5*DZ./(DX.*DY)./kz;			% same for z direction

Cx=1./(Rx(:,1:end-1,:)+Rx(:,2:end,:));		% conductances x-direction
Cy=1./(Ry(1:end-1,:,:)+Ry(2:end,:,:));		% same for y direction
Cz=1./(Rz(:,:,1:end-1)+Rz(:,:,2:end));		% same for z direction

A=sparse([Ix(:);Jx(:);Iy(:);Jy(:);Iz(:);Jz(:)],...
         [Jx(:);Ix(:);Jy(:);Iy(:);Jz(:);Iz(:)],...
        -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],Ncel,Ncel,7*Ncel);	% matrix without main diag
Adiag= -sum(A,2);	% diagonal, kept saparately

% Boundary conditions, just Q and Fixed Heads right now
isFxHd=~isnan(FH); FH(~isFxHd)=0;

PCG=1;
if PCG
	AA=spdiags(isFxHd*HUGE + ~isFxHd.*Adiag,0,A);						% Total matrix to solve
	M=cholinc(AA,1e-3);															% cholesky precoditioner
	Phi=pcg(AA,(isFxHd.*FH*HUGE + ~isFxHd.*Q),1e-25,100,M',M);		% solve the system
else
	Phi=spdiags(isFxHd*HUGE + ~isFxHd.*Adiag,0,A)\(isFxHd.*FH*HUGE + ~isFxHd.*Q);	% may cost less time
end

Q=spdiags(Adiag,0,A)*Phi;

Q  =reshape(Q,Ny,Nx,Nz);
Phi=reshape(Phi,Ny,Nx,Nz);

if nargout>2				% flow across cell faces
	QRF=-diff(Phi,1,2).*Cx;
	QFF=-diff(Phi,1,1).*Cy;
	QLF=-diff(Phi,1,3).*Cz;
end
