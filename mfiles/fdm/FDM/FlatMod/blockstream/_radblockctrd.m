function [Phi,Q,Psiy,Psix]=radblockctrd(x,y,kx,ky,FH,FQ)
% function [Phi,Q,Psiy,Psix]=radblockctrd(x,y,kx,ky,FH,FQ)
% Defines and solves a 2D radial block-centred finite difference model
% x,y mesh coordinates, kx,ky conductivities, FH=fixed heads (NaN for ordinary points) FQ=fixed nodal flows
% Phi is computed heads, Q computed nodal balances
% Psiy is stream function with vertical branch cuts
% Psix is stream function with horizontal branch cuts
% Use thin outer cells to get true boundary flow.
% Is using given heads in outer cells, only half the width will be dealt with in the computation.
% TO 991017, TO 000525/000530

HUGE=1e20;			% Een heel groot getal (zie verderop))))
TINY=1e-20;			% Een heel klein getal (ter verkoming van deling door nul, zie verderop)

x=x(:)';  xM=(x(1:end-1)+x(2:end))/2;			% celmiddens
y=y(:);   yM=(y(1:end-1)+y(2:end))/2;			% celmiddens

Nx=length(xM);  Ny=length(yM);					% cellen tellen

dx=diff(x);        dx(find(dx==0))=TINY;		% celafmetingen
dy=abs(diff(y));   dy(find(dy==0))=TINY;		% celafmetingen

if isempty(FH),  FH=NaN*zeros(Ny,Nx); end; FH=FH(:);		% Eventueel matrix voor gegeven stijghoogten
if isempty(FQ),  FQ=    zeros(Ny,Nx); end; FQ=FQ(:);		% Eventueel matrix voor gegeven volumestromen

% node numbering
Nodes = reshape([1:Nx*Ny],Ny,Nx);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

[DX,DY]=meshgrid(dx,dy);						% Volledige matrix aanmaken
x =ones(Ny,1)*x;
xM=ones(Ny,1)*xM;

% resistance to flow in cells, radial
RXR=(log(xM./x(:,1:end-1))+TINY)./(2*pi*kx.*DY); RXR=RXR(:,2:end  );		% Weerstand radiaal
RXL=(log(x(:,2:end  )./xM)+TINY)./(2*pi*kx.*DY); RXL=RXL(:,1:end-1);		% Weerstand radiaal
ex=1./(RXL+RXR);									% Conductance horizontaal

% resistance vertical, in rings
Aring=pi*(x(:,2:end  ).^2-x(:,1:end-1).^2); Aring(find(x(:,2:end)==x(:,1:end-1)))=TINY;
RY=0.5*DY./(Aring.*ky);							% Weerstand verticaal
RYT=RY(1:end-1,:);
RYB=RY(2:end  ,:);
ey=1./(RYT+RYB);									% Conductance verticaal

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [ex(:);ex(:);ey(:);ey(:)],Ny*Nx,Ny*Nx,5*Ny*Nx);	% Systeemmatrix
Adiag= -sum(A,2);									% Diagonaal van de systeemmatrix

% Boundary conditions, just FQ and Fixed Heads right now
isFxHd=HUGE* ~isnan(FH);
FH(~isFxHd)=0;

Phi=spdiags(isFxHd + ~isFxHd.*Adiag,0,A)\(isFxHd.*FH + ~isFxHd.*FQ);	% Berekening stijghoogten

Q=spdiags(Adiag,0,A)*Phi;						% Berekening volumestromen (injecties in de cellen)

Q  =reshape(Q,Ny,Nx);							% Terugvouwen in de vorm van het oorspronkelijke netwerk
Phi=reshape(Phi,Ny,Nx);							% Terugvouwen in de vorm van het oorspronkelijke netwerk

qx= ex.*diff(Phi,1,2); Qx=qx(2:end-1,:);	% volumestroom over de verticale celwanden
qy= ey.*diff(Phi,1,1); Qy=qy(:,2:end-1);	% volumestroom over de horizontale celwanden

Psix=zeros(Ny+1,Nx+1);							% stroomfunctie met horizontale branchcuts
Psix(1:end-1,2:end-1)=Psix(1:end-1,2:end-1)+flipud(cumsum(flipud(qx),1));

Psiy=zeros(Ny+1,Nx+1);							% stroomfunctie met verticale branchcuts
Psiy(2:end-1,2:end)=Psiy(2:end-1,2:end)+cumsum(qy,2);