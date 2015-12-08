function [Phi,Q,Psi]=blkstrm(x,y,kx,ky,FH,Q)
% function [Phi,Q,Psi]=blkctrd(x,y,kx,ky,FH,Q)
% 2D block-centred steady state finite difference model with stream function
% x,y mesh coordinates, kx,ky conductivities, FH=fixed heads (NaN for ordinary points) Q=fixed nodal flows
% Phi,Q computed heads and cell balances
% Psi, stream function with vertical branch cuts at the corners of the cells first and last columns zero
% and should be omitted when contouring, stream function at bottom of model=0
% TO 991017  TO 000530 001026

HUGE=1e+20;
TINY=1e-10;

x=x(:)'; Nx=length(x)-1; dx=abs(diff(x));	dx=max(TINY,dx);
y=y(:);  Ny=length(y)-1; dy=abs(diff(y));	dy=max(TINY,dy);

if isempty(FH), error('Need a fixed heads array, use help blkctrd'); end
if isempty( Q),   Q=    zeros(Ny,Nx); end;  Q= Q(:);	% Generate fixed flow matrix if necessary

% node numbering
Nodes = reshape([1:Nx*Ny],Ny,Nx);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

[DX,DY]=meshgrid(dx,dy);						% Generate complete mesh

RX=0.5*DX./DY./kx; ex=1./(RX(:,1:end-1)+RX(:,2:end));	% Resistance and conductance
RY=0.5*DY./DX./ky; ey=1./(RY(1:end-1,:)+RY(2:end,:)); % Resistance and conductance

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [ex(:);ex(:);ey(:);ey(:)],Ny*Nx,Ny*Nx,5*Ny*Nx);		% system matrix
Adiag= -sum(A,2);																% Main diagonal of system matrix

isFxHd=HUGE* ~isnan(FH(:));	FH(~isFxHd)=0; 							% head boundary condition

Phi=spdiags(isFxHd + ~isFxHd.*Adiag,0,A)\(isFxHd.*FH(:) + ~isFxHd.*Q);		% Compute heads

Phi=reshape(Phi,Ny,Nx);							% Always reshape Phi back into shape of original model
if nargout>1
	Q=spdiags(Adiag,0,A)*Phi(:);						% Cell flows only if wanted
   Q  =reshape(Q,Ny,Nx);							% reshape Q back into model shape
	if nargout>2,											% vertical branch cuts
		qx=ex.*diff(Phi,1,2);							% Flow across vertical cell walls
		Psi=zeros(Ny+1,Nx+1);							% Stream function matrix with horizontal branch cuts
		Psi(1:end-1,2:end-1)=Psi(1:end-1,2:end-1)+flipud(cumsum(flipud(qx),1));
   end
end

