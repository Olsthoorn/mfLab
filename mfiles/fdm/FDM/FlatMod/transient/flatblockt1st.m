function   [Phi,Q,QRF,QFF]=flatblockt1st(x,y,dt,kx,ky,S,FH,FQ,Cghb,IH,theta)
% function [Phi,Q,QRF,QFF]=flatblockt1st(x,y,t,kx,ky,S,FH,FQ,IH,[Cghb,theta])
% Single time step transient single-layer block-centered groundwater finite difference model
% x,y,t coordinates,
% kx,ky conductivities, S=storage coef., use NaNs in one of these or IH to make cells inactive
% FH=fixed heads (NaN for ordinary points)
% FQ=fixed nodal flows
% IH=initial heads
% Cghb is general head conductance applied to the FH.
% Cghb can be used as ibound.
%		If ~isnan(FH) & ~isnan(Cghb) 				--> FH and Cghb are applied,
%		IF ~isnan(FH) &  isnan(Cghb) |Chgb==0 	--> FH treated as true fixed head, if <0  
% Phi is computed heads, Q computed nodal balances
% QRF=Q right face (between cell columns)
% QFF=Q front face (between cell rows)
% TO 991017, 010809 010829, 010904

%persistent Nx Ny A Adiag IA								% use this if you want to skip maxtrix computation after call 1

HUGE=1e20;

if nargin< 9 | isempty(Cghb),   Cghb=NaN*FH; end
if nargin<10 | isempty(theta),  theta=0.67;  end	% impliciteness

x=x(:)';  xM=(x(1:end-1)+x(2:end))/2;
y=y(:);   yM=(y(1:end-1)+y(2:end))/2;

Nx=length(xM);  Ny=length(yM);
dx=diff(x);    dy=diff(y); if y(end)<y(1), dy=-dy; end

S=S.*(dy*dx);

% Fixed Heads and General Head Boundaries
Cghb(find(Cghb<=0 | isnan(Cghb)))=   0;		% illegal GHConductance, becomes fixed head
Cghb(find(~isnan(FH) & Cghb==0)) =HUGE;	% truly fixed head nodes
FH  (find(isnan(FH  )))=0;								% we now have the correct conductance everywhere, so NaNs may go

%if isempty(A),											% use this to skip computing matrix after run one
%fprintf('make matrix\n');

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

IA=find(~isnan(S(:)+kx(:)+ky(:)+IH(:)));			% IA are the active cells, kx, ky, IH, S may have NaN to indicate inactive cells
Adiag= -sum(A(IA,IA),2);								% Sum over the active cells only

%end															% use this to skip matrix computation after run 1 
%fprintf('.');

Phi=NaN*zeros(size(FH(:)));
Q  =NaN*zeros(size(FH(:)));

phi =spdiags(Cghb(IA)+S(IA)./(dt*theta)+Adiag,0,A(IA,IA))\(Cghb(IA).*FH(IA) + FQ(IA) + S(IA)./(dt*theta).*IH(IA));

Phi(IA)=IH(IA)+(phi-IH(IA))/theta;
Q  (IA)=spdiags(Adiag,0,A(IA,IA))*Phi(IA);   

Q  =reshape(Q,Ny,Nx);
Phi=reshape(Phi,Ny,Nx);

if nargout>2,									% compute right face and front face flows
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
     