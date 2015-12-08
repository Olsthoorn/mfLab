function   [Phi,Q,QRF,QFF]=fd2bt(x,y,t,kx,ky,S,FH,FQ,IH,Cghb,theta)
% FUNCTION [Phi,Q,QRF,QFF]=fd2bt(x,y,t,kx,ky,S,FH,FQ,IH,[Cghb,theta])
% Transient explicit single layer groundwater finite difference model
% Defines and solves a 2D block-centred finite difference model
% x,y,t coordinates,  x must be a row vector, y must be a column vector or a 2d matrix with depth vertical (first dimension in matlab)
% kx,ky conductivities, S=specific storage coef.,
% FH=fixed heads (NaN for ordinary points)
% FQ=fixed nodal flows
% IH=initial heads
% Cghb is general head conductance applied to the FH.
% Cghb can be used as ibound.
%		If ~isnan(FH) & ~isnan(Cghb) 				--> FH and Cghb are applied,
%		IF ~isnan(FH) &  isnan(Cghb) |Chgb==0 	--> FH treated as true fixed head, if <0  
% Phi is computed heads, Q computed nodal balances
% QRF=Q right face, including a zeros column to the right, like in MODFLOW
% QFF=Q front face, inclufin a zeros column to the front, like in MODFLOW
% TO 991017, 010809 010829

HUGE=1e20;

t=t(:)';  dt=diff(t);
xM=(x(1:end-1  )+x(2:end  ))/2;
yM=(y(1:end-1,:)+y(2:end,:))/2;

Nx=size(xM,2);  Ny=size(yM,1); Nt=length(t);
dx=abs(diff(x,1,2));
dy=abs(diff(y,1,1));

if isempty(FH),  FH=NaN*zeros(Ny,Nx); end; FH=FH(:);
if isempty(FQ),  FQ=    zeros(Ny,Nx); end; FQ=FQ(:);
if nargin< 9 | isempty(Cghb),   Cghb=NaN*FH; end
if nargin<10 | isempty(theta),  theta=0.67;  end	% impliciteness

if size(dy,2)==size(dx,2),              % full matrix for dy but row matrix for dx (always dx must be a row matrix)
    DX=ones(size(dy(:,1)))*dx;
    DY=dy;
else
    [DX,DY]=meshgrid(dx,dy);
end

S=S.*DY.*DX; S=S(:);

% Fixed Heads and General Head Boundaries
Cghb(find(Cghb<=0 | isnan(Cghb)))	=   0;		% illegal GHConductance, becomes fixed head
Cghb(find(~isnan(FH(:)) & Cghb(:)==0)) =HUGE;	% truly fixed head nodes
FH  (find(isnan(FH  )))=0;								% we now have the correct conductance everywhere, so NaNs may go

% node numbering
Nodes = reshape([1:Nx*Ny],Ny,Nx);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

RH=0.5*DX./DY./kx;  ex=1./(RH(:,1:end-1)+RH(:,2:end));
RV=0.5*DY./DX./ky;  ey=1./(RV(1:end-1,:)+RV(2:end,:));

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [ex(:);ex(:);ey(:);ey(:)],Ny*Nx,Ny*Nx,5*Ny*Nx);
Adiag= -sum(A,2);


Phi=NaN*zeros(length(FH(:)),length(t));
Q  =NaN*zeros(length(FH(:)),length(t));

%Phi(:,1)=spdiags(Cghb + Adiag,0,A)\(Cghb.*FH + FQ(:,1));		% steady state solution
Phi(:,1)=IH(:);
Q  (:,1)=spdiags(Adiag,0,A)*Phi(:,1);

% dealing with inactive cells, selecting only the active cells
IA=find(~isnan(S(:)+kx(:)+ky(:)+IH(:)));		% Inactive cells are earmarked by NaNs in kx,ky or S of IH
A=A(IA,IA); Adiag= -sum(A,2); FH=FH(IA); FQ=FQ(IA); Cghb=Cghb(IA); S=S(IA);	IH=IH(IA);	% only active cells

tic
fprintf('time steps:');
for it=1:length(dt)
   if ~rem(it,50), fprintf('\n'); end
   phi=spdiags(Cghb + S./(dt(it)*theta) + Adiag,0,A)\(Cghb.*FH + FQ + S./(dt(it)*theta).*Phi(IA,it));
   fprintf('.');
   Phi(IA,it+1)=Phi(IA,it)+(phi-Phi(IA,it))/theta;
   Q  (IA,it+1)=spdiags(Adiag,0,A)*Phi(IA,it);   
end
fprintf('ready in %f seconds\n',toc);

Q  =reshape(Q,Ny,Nx,Nt);
Phi=reshape(Phi,Ny,Nx,Nt);

if nargout>2
    signx=sign(x(1,end)-x(1,1));
    signy=sign(y(end,1)-y(1,1));							% regular direction, y axis downward if rows numbered upard
    QRF=-signx*diff(Phi,1,2);
    QFF=-signy*diff(Phi,1,1);
    for iz=1:size(Phi,3)
        QRF(:,:,iz)=QRF(:,:,iz).*ex;
        QFF(:,:,iz)=QFF(:,:,iz).*ey;
    end
end