function [Phi,Qt,Qx,Qy,Qz,Qs]=fdm3t(x,y,z,t,kx,ky,kz,S,IH,IBOUND,FQ)
%FDM3T a 3D block-centred transient finite difference model
%
% Exemple:
%    [Phi,Q,Qx,Qy,Qz,Qs]=fdm3t(x,y,z,t,kx,ky,kz,S,IH,IBOUND,FQ)
%
% INPUT:
%  x(NX+1)         x-coordinate of mesh/grid
%  y(NY+1)         y-coordinate of mesh/grid
%  z(NZ+1)         z-coordinate of mesh/grid
%  t(Nt+1)         time for output. 0 will be added if necessary
%  kx(NY,NX,NZ)    conductivity in x-direction
%  ky(NY,NX,NZ)    same in y direction
%  kz(NY,NX,NZ)    same in z direction
%  S(NY,NX,NZ)     primary storage (S+Sy)
%  IH(NY,NX,NZ)    initial head
%  IBOUND(NY,NX,NZ) determines fixed head (-1) and inactive cell (0) lccations,
%  Q=fixed nodal flows
%
% OUTPUT
%  Phi(NY,NX,NZ,Nt+1) computed heads with Phi(NY,NX,1) initial heads for t=0
%, Qt(NY,NX,  NZ,Nt)    computed total cell balance during timestep it
%  Qx(NY,NX-1,NZ,Nt)  hor.  cell face flow in x-direction positive along increasing col number
%  Qy(NY-1,NX,NZ,Nt)  vert. cell face flow in y-directin  positive along increasing row number
%  Qs(NY,NX,NZ-1,Nt)    storage change of node during timestep it
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301 100307

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

theta=0.67;  % implicitness

if isempty(ky), ky=kx; end

if t(1)~=0, t=unique([0 abs(t)]); else t=unique(abs(t)); end     % Initial time is always zero

sx=sign(x(end)-x(1)); x=unique(x); dx=diff(x); NX=length(dx); %xm=0.5*(x(1:end-1)+x(2:end));
sy=sign(y(end)-y(1)); y=unique(y); dy=diff(y); NY=length(dy); %ym=0.5*(y(1:end-1)+y(2:end));
sz=sign(z(end)-z(1)); z=unique(z); dz=diff(z); NZ=length(dz); %zm=0.5*(z(1:end-1)+z(2:end));

if sx<0
    kx=kx(:,end:-1:1,:);
    ky=ky(:,end:-1:1,:);
    kz=kz(:,end:-1:1,:);
    S =S( :,end:-1:1,:);
    IH=IH(:,end:-1:1,:);
    IBOUND=IBOUND(:,end:-1:1,:);
    FQ=FQ(:,end:-1:1,:);
end
if sy<0
    kx=kx(end:-1:1,:,:);
    ky=ky(end:-1:1,:,:);
    kz=kz(end:-1:1,:,:);
    S =S( end:-1:1,:,:);
    IH=IH(end:-1:1,:,:);
    IBOUND=IBOUND(end:-1:1,:,:);
    FQ=FQ(end:-1:1,:,:);
end
if sz<0
    kx=kx(:,:,end:-1:1);
    ky=ky(:,:,end:-1:1);
    kz=kz(:,:,end:-1:1);
    S =S( :,:,end:-1:1);
    IH=IH(:,:,end:-1:1);
    IBOUND=IBOUND(:,:,end:-1:1);
    FQ=FQ(:,:,end:-1:1);
end

[DX,DY,DZ]=meshgrid(dx,dy,dz); NODES=NY*NX*NZ;

t=unique([0; t(:)]);  dt=    diff(t);  Nt=length(dt);

Nodes = reshape(1:NODES,NY,NX,NZ);               % Node numbering
IE=Nodes(:,2:end,:);   IW=Nodes(:,1:end-1,:);
IS=Nodes(2:end,:,:);   IN=Nodes(1:end-1,:,:);
IT=Nodes(:,:,2:end);   IB=Nodes(:,:,1:end-1);

% resistances and conducctances
RX=0.5*(DX./(DY.*DZ))./kx; Cx=1./(RX(:,1:end-1,:)+RX(:,2:end,:));
RY=0.5*(DY./(DZ.*DX))./ky; Cy=1./(RY(1:end-1,:,:)+RY(2:end,:,:));
RZ=0.5*(DZ./(DX.*DY))./kz; Cz=1./(RZ(:,:,1:end-1)+RZ(:,:,2:end));

Cs=DX.*DY.*DZ.*S;

Cs=Cs(:);  % storage conductacne when devided by dt*theta

A=sparse([IE(:);IW(:);IN(:);IS(:);IT(:);IB(:)],...
         [IW(:);IE(:);IS(:);IN(:);IB(:);IT(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],...
         NODES,NODES,7*NODES);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct =Nodes((kx>0 | ky>0 | kz>0 | S>0) &  IBOUND  );              % active cells
I    =Nodes((kx>0 | ky>0 | kz>0 | S>0) &  IBOUND>0); % active cells but not fixed heads = cells with heads to be computed
Ifh  =Nodes((kx>0 | ky>0 | kz>0 | S>0) &  IBOUND<0); % active cells with fixed heads

Phi=NaN(NODES,Nt+1);  % allocate space to store the entire head matrix
Qt =NaN(NODES,Nt);    % allocate memory for Qt
Qs =NaN(NODES,Nt);    % allocate memory for Qs
Qx =NaN(NY,NX-1,NZ,Nt);  % allocate memory for Qx
Qy =NaN(NY-1,NX,NZ,Nt);  % allocate memory for Qy
Qz =NaN(NY,NX,NZ-1,Nt);  % allocate memory for Qy


Phi(IAct,1)=IH(IAct);  % store initial head at Phi(:,:,1)
FQ=FQ(:);
Fi =IH(:);            % head computed in this time step
if any(I(:))
    fprintf('Time steps:');
    for it=1:Nt
        fprintf('.');
        Fi (I)         =spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I)-A(I,Ifh)*IH(Ifh)+Cs(I).*Phi(I,it)/dt(it)/theta); % solve
        Phi(IAct ,it+1)=Fi(IAct)/theta-(1-theta)/theta*Phi(IAct,it);
        Qt (IAct ,it  )=spdiags(Adiag(IAct),0,A(IAct,IAct))*Fi(IAct);
        Qs (IAct ,it  )=-Cs(IAct).*(Phi(IAct,it+1)-Phi(IAct,it))/dt(it);   % Storage in time step m3 for cell
        Qx (:,:,:,it  )=-Cx.*diff(reshape(Fi,NY,NX,NZ),1,2);	 % Flow across horizontal cell faces m3/d per cell
        Qy (:,:,:,it  )=-Cy.*diff(reshape(Fi,NY,NX,NZ),1,1);   % Flow across vertical cell faces, m3/d per cell
        Qz (:,:,:,it  )=-Cz.*diff(reshape(Fi,NY,NX,NZ),1,3);   % Flow across vertical cell faces, m3/d per cell
    end
    fprintf('%d time steps done\n',Nt);
end
Phi=reshape(Phi,NY,NX,NZ,Nt+1);                   % NaN if inactive
Qt =reshape(Qt ,NY,NX,NZ,Nt  ); Qt(isnan(Qt))=0;  % 0 if inactive
Qs =reshape(Qs ,NY,NX,NZ,Nt  ); Qs(isnan(Qs))=0;  % 0 if inactive
Qx(isnan(Qx))=0; Qy(isnan(Qy))=0; Qz(isnan(Qz))=0; % 0 if inactive

if sx<0
    Phi=Phi(:,end:-1:1,:,:);
    Qt =Qt( :,end:-1:1,:,:);
    Qs =Qt( :,end:-1:1,:,:);
    Qx =Qx( :,end:-1:1,:,:);
    Qy =Qy( :,end:-1:1,:,:);
    Qz =Qz( :,end:-1:1,:,:);
end
if sy<0
    Phi=Phi(end:-1:1,:,:,:);
    Qt =Qt( end:-1:1,:,:,:);
    Qs =Qt( end:-1:1,:,:,:);
    Qx =Qx( end:-1:1,:,:,:);
    Qy =Qy( end:-1:1,:,:,:);
    Qz =Qz( end:-1:1,:,:,:);
end
if sz<0
    Phi=Phi(:,:,end:-1:1,:);
    Qt =Qt( :,:,end:-1:1,:);
    Qs =Qt( :,:,end:-1:1,:);
    Qx =Qx( :,:,end:-1:1,:);
    Qy =Qy( :,:,end:-1:1,:);
    Qz =Qz( :,:,end:-1:1,:);
end
Phi=Phi(:,:,:,2:end);
