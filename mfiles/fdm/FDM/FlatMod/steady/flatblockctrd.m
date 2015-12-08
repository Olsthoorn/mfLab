function [Phi,Q,QRF,QFF]=flatblockctrd(x,y,kx,ky,FH,Q,Cghb)
% function [Phi,Q,QRF,QFF]=flatblockctrd(x,y,kx,ky,FH,Q,Cghb)
% Defines and solves a 2D block-centred finite difference model
% x,y mesh coordinates, kx,ky conductivities, FH=fixed heads (NaN for ordinary points) Q=fixed nodal flows
% Cghb is general head conductance applied to the FH.
% FH~=NaN, then if Cghb~=NaN, apply Cghb as conductance, else apply fixed head (Cghb becomes HUGE)
% Phi is computed heads, Q computed nodal balances
% QRF=Q right face, including a zeros column to the right, like in MODFLOW
% QFF=Q front face, inclufin a zeros column to the front, like in MODFLOW
% TO 991017, 010809

HUGE=1e20;

Eghb=~isnan(FH(:))*HUGE; if nargin>6 & ~isempty(Cghb), I=find(~isnan(Cghb)); Eghb(I)=Cghb(I); end

x=x(:)';  xM=(x(1:end-1)+x(2:end))/2;
y=y(:);   yM=(y(1:end-1)+y(2:end))/2;

Nx=length(xM);  Ny=length(yM);

dx=diff(x);    dy=diff(y); if y(end)<y(1), dy=-dy; end

if isempty(FH),  FH=NaN*zeros(Ny,Nx); end; FH=FH(:);
if isempty( Q),   Q=    zeros(Ny,Nx); end;  Q= Q(:);

% node numbering
Nodes = reshape([1:Nx*Ny],Ny,Nx);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

[DX,DY]=meshgrid(dx,dy);

RH=0.5*DX./DY./kx;  ex=1./(RH(:,1:end-1)+RH(:,2:end));
RV=0.5*DY./DX./ky;  ey=1./(RV(1:end-1,:)+RV(2:end,:));

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [ex(:);ex(:);ey(:);ey(:)],Ny*Nx,Ny*Nx,5*Ny*Nx);
Adiag= -sum(A,2);

% Boundary conditions, just Q and Fixed Heads right now
Eghb(find(isnan(Eghb)))=0;
FH  (find(isnan(FH  )))=0;

Phi=spdiags(Eghb + Adiag,0,A)\(Eghb.*FH + Q);

Q=spdiags(Adiag,0,A)*Phi;

Q  =reshape(Q,Ny,Nx);
Phi=reshape(Phi,Ny,Nx);

if nargout>2
   if x(end)<x(1)
      QRF=+diff(Phi,1,2).*ex;
   else											% regualar direction, x axis upward to right
      QRF=-diff(Phi,1,2).*ex;
   end
   
   if y(end)<y(1),							% regular direction, y axis downward if rows numbered upard
      QFF=+diff(Phi,1,1).*ey;
   else
      QFF=-diff(Phi,1,1).*ey;
   end
end
      