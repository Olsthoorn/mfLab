function   [Phi,Q,QRF,QFF]=flatblockctrdt(x,y,t,kx,ky,S,FH,FQ,Cghb,IH,theta,solver)
% function [Phi,Q,QRF,QFF]=flatblockctrdt(x,y,t,kx,ky,S,FH,FQ,IH,[Cghb,theta,solver])
% Transient explicit single layer groundwater finite difference model
% Defines and solves a 2D block-centred finite difference model
% x,y,t coordinates,
% kx,ky conductivities, S=storage coef.,
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
x=x(:)';  xM=(x(1:end-1)+x(2:end))/2;
y=y(:);   yM=(y(1:end-1)+y(2:end))/2;

Nx=length(xM);  Ny=length(yM); Nt=length(t);
dx=diff(x);    dy=diff(y); if y(end)<y(1), dy=-dy; end

if isempty(FH),  FH=NaN*zeros(Ny,Nx); end; FH=FH(:);
if isempty(FQ),  FQ=    zeros(Ny,Nx); end; FQ=FQ(:);
if nargin< 9 | isempty(Cghb),   Cghb=NaN*FH; end
if nargin<10 | isempty(theta),  theta=0.67;  end	% impliciteness
if nargin<11 | isempty(solver), solver=1;    end	% solver 1 = direct, 2 = pcg with choleski preconditioner

S=S.*(dy*dx);  S=S(:);

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

[DX,DY]=meshgrid(dx,dy);

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
pcgflag=1; TOL=1e-6; MAXIT=30;
for it=1:length(dt)
   if ~rem(it,50), fprintf('\n'); end
   switch solver
      case 1,			% direct solution
         phi		  =spdiags(Cghb + S./(dt(it)*theta) + Adiag,0,A)\(Cghb.*FH + FQ + S./(dt(it)*theta).*Phi(IA,it));
  	   	fprintf('.');
      case 2,			% pcg with precoditioner
			cnt=0;
			while 1
      		if pcgflag,	% preconditioner computation only if necessary
         		R = cholinc(spdiags(Cghb + S./(dt(it)*theta) + Adiag,0,A),1e-3);
	         	fprintf('R');
	   	   end
				[phi,pcgflag] = pcg(spdiags(Cghb + S./(dt(it)*theta) + Adiag,0,A),(Cghb.*FH + FQ + S./(dt(it)*theta).*Phi(IA,it)),TOL,MAXIT,R',R);
      	   fprintf('%d',pcgflag);
	      	if ~pcgflag,
   	      	cnt=0;
	      	   break;
		      else, 	% if not converging after preconditioning
   		      cnt=cnt+1;
      		   if cnt>1, break; end
            end
         end
      otherwise,			% illegal solver
         error('use 1 or 2 for the solver');
      end
   Phi(IA,it+1)=Phi(IA,it)+(phi-Phi(IA,it))/theta;
   Q  (IA,it+1)=spdiags(Adiag,0,A)*Phi(IA,it);   
end
fprintf('ready in %f seconds\n',toc);

Q  =reshape(Q,Ny,Nx,Nt);
Phi=reshape(Phi,Ny,Nx,Nt);

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
      