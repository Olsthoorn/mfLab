function [Phi,Q]=FD3DMC(x,y,z,kx,ky,kz,FH,Q)
% function [Phi,Q]=FD3DMC(x,y,z,kx,ky,kz,FH,Q)
% Defines and solves a 3D mesh-centred finite difference model
% x,y,z mesh coordinates, kx,ky,kz conductivities, FH=fixed heads (NaN for ordinary points) Q=fixed nodal flows
% Phi is computed heads, Q computed nodal balances
%
% TO 991119, TO 000120 TO 010520

HUGE=1e20;

x=x(:)';	Nx=length(x);	dx=abs(diff(x));
y=y(:);	Ny=length(y);	dy=abs(diff(y)); 
z=z(:)';	Nz=length(z);	dz=abs(diff(z));	Nod=Nx*Ny*Nz; [DX,DY,DZ]=meshgrid(dx,dy,dz);

if isempty(FH),  FH=NaN*zeros(Ny,Nx,Nz); end; FH=FH(:);
if isempty( Q),   Q=    zeros(Ny,Nx,Nz); end;  Q= Q(:);

% node numbering
Nodes = reshape([1:Nx*Ny*Nz],Ny,Nx,Nz);
Ix=Nodes(:,1:end-1,:); Jx=Nodes(:,2:end  ,:);
Iy=Nodes(1:end-1,:,:); Jy=Nodes(2:end  ,:,:);
Iz=Nodes(:,:,1:end-1); Jz=Nodes(:,:,2:end  );

Ex=0.25*kx.*DY.*DZ./DX;
Ey=0.25*ky.*DX.*DZ./DY;
Ez=0.25*kz.*DX.*DY./DZ;

Cx=zeros(Ny,Nx-1,Nz);
Cx(1:end-1,:,1:end-1)=Cx(1:end-1,:,1:end-1)+Ex;
Cx(2:end  ,:,1:end-1)=Cx(2:end  ,:,1:end-1)+Ex;
Cx(1:end-1,:,2:end  )=Cx(1:end-1,:,2:end  )+Ex;
Cx(2:end  ,:,2:end  )=Cx(2:end  ,:,2:end  )+Ex;

Cy=zeros(Ny-1,Nx,Nz);
Cy(:,1:end-1,1:end-1)=Cy(:,1:end-1,1:end-1)+Ey;
Cy(:,2:end  ,1:end-1)=Cy(:,2:end  ,1:end-1)+Ey;
Cy(:,1:end-1,2:end  )=Cy(:,1:end-1,2:end  )+Ey;
Cy(:,2:end  ,2:end  )=Cy(:,2:end  ,2:end  )+Ey;

Cz=zeros(Ny,Nx,Nz-1);
Cz(1:end-1,1:end-1,:)=Cz(1:end-1,1:end-1,:)+Ez;
Cz(1:end-1,2:end  ,:)=Cz(1:end-1,2:end  ,:)+Ez;
Cz(2:end  ,1:end-1,:)=Cz(2:end  ,1:end-1,:)+Ez;
Cz(2:end  ,2:end  ,:)=Cz(2:end  ,2:end  ,:)+Ez;

A=sparse([Ix(:);Jx(:);Iy(:);Jy(:);Iz(:);Jz(:)],...
         [Jx(:);Ix(:);Jy(:);Iy(:);Jz(:);Iz(:)],...
        -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],Nod,Nod,7*Nod);
Adiag= -sum(A,2);

% Boundary conditions, just Q and Fixed Heads right now
isFxHd= ~isnan(FH); FH(~isFxHd)=0;

PCG=0;
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
