function [Phi,Q,Psiy,Psix]=blkstrm(x,y,kx,ky,FH,Q)
% function [Phi,Q,Psiy,Psix]=blkctrd(x,y,kx,ky,FH,Q)
% 2D block-centred steady state finite difference model
% x,y mesh coordinates, kx,ky conductivities, FH=fixed heads (NaN for ordinary points) Q=fixed nodal flows
% Phi,Q computed heads and cell balances
% Psiy, stream function with vertical branch cuts		(at the corners of the cells, first and last column=0,
% 		stream function at bottom of model=0)
% Psix, stream function with horizontal branch cuts	(at the corners of the cells, first and last row=0,
%		stream function at left column=0)
% TO 991017  TO 000530 001026

HUGE=1e20;

x=x(:)'; Nx=length(x)-1; dx=diff(x);
y=y(:);  Ny=length(y)-1; dy=abs(diff(y));

if isempty(FH), error('Need a fixed heads array, use help blkctrd'); end
if isempty( Q),   Q=    zeros(Ny,Nx); end;  Q= Q(:);	% Generate fixed flow matrix if necessary

% node numbering
Nodes = reshape([1:Nx*Ny],Ny,Nx);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

[DX,DY]=meshgrid(dx,dy);						% Generate complete mesh

warning off
RX=0.5*DX./DY./kx; ex=1./(RX(:,1:end-1)+RX(:,2:end));	% Resistance and conductance
RY=0.5*DY./DX./ky; ey=1./(RY(1:end-1,:)+RY(2:end,:));	% Resistance and conductance
warning on

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [ex(:);ex(:);ey(:);ey(:)],Ny*Nx,Ny*Nx,5*Ny*Nx);		% system matrix
Adiag= -sum(A,2);																% Main diagonal of system matrix

isFxHd=HUGE* ~isnan(FH(:));	FH(~isFxHd)=0; 						% head boundary condition

Phi=spdiags(isFxHd + ~isFxHd.*Adiag,0,A)\(isFxHd.*FH(:) + ~isFxHd.*Q);		% Compute heads

Phi=reshape(Phi,Ny,Nx);							% Always reshape Phi back into shape of original model
if nargout>1
	Q=spdiags(Adiag,0,A)*Phi(:);					% Cell flows only if wanted
   Q  =reshape(Q,Ny,Nx);							% reshape Q back into model shape
	if nargout>2,											% vertical branch cuts
		qy=ey.*diff(Phi,1,1);							% Flow over de horizontale celwanden
		Psiy=zeros(Ny+1,Nx+1);							% Idem verticale branchcuts
		Psiy(2:end-1,2:end)=Psiy(2:end-1,2:end)+cumsum(qy,2);
		if nargout>3,											% horizontal branch cuts
			qx=ex.*diff(Phi,1,2);							% Flow across vertical cell walls
			Psix=zeros(Ny+1,Nx+1);							% Stream function matrix with horizontal branch cuts
			Psix(1:end-1,2:end-1)=Psix(1:end-1,2:end-1)+flipud(cumsum(flipud(qx),1));
   	end
   end
end

