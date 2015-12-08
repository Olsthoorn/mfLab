function [Phi,Qt,Qx,Qy,Qs]=fdm2tNC(xGr,yGr,t,Kx,Ky,S,IBOUND,IH,FQ,mode)
%FDM2TNC a 2D node-centred transient finite difference model
%
% Example:
%    [Phi,Q,Qx,Qy,Qs]=fdm2ctNodeCtrd(xGr,yGr,t,Kx,c,Ky,S,IH,FH,FQ [,radial])
%
% INPUT:
%  xGr(Nx+1)       xGr-coordinate of mesh/grid
%  yGr(Ny+1)       yGr-coordinate of mesh/grid
%  t(Nt+1)         time for output. 0 will be added if necessary
%  Kx(Ny,Nx)       conductivity in x-direction
%  Ky(Ny,Nx)       same in y direction, Ky=Kx if Ky==[]
%  S(Ny,Nx)        primary storage (S+Sy)
%  IH(Ny,Nx)       initial head
%  FH(Ny,Nx)       fixed heads (NaN for ordinary points), Q=fixed nodal flows
%  Radial          Arbitrary input caused model to assume axial flow with
%                 xGr=r and FQ are ring flows while als flows have dimension
%                 L3/T instead of L2/t. For clearness use 'R' or 'Radial'
%                 as input at the function call.
% OUTPUT
%  Phi(Ny,Nx,Nt+1) computed heads with Phi(Ny,Nx,1) initial heads for t=0
%, Qt(Ny,Nx,Nt)    computed total cell balance during timestep it
%  Qx(Ny,Nx-1,Nt)  hor.  cell face flow in x-direction positive along increasing col number
%  Qy(Ny-1,Nx,Nt)  vert. cell face flow in y-direction positive along increasing row number
%  Qs(Ny,Nx,Nt)    storage change of node during timestep it
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301 120111

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

theta= 0.67; % implicitness

[xGr,~,xm,~,dx,dy,Nx,Ny]=modelsize(xGr,yGr);
    
t=unique(t(:));   dt= diff(t);   Nt=length(dt);

if all(IBOUND==0), error('IBOUND all zeros (inactive), nothing to compute!'); end
if size(Kx,2)==1, Kx=Kx*ones(size(xm)); end
if size(Ky,2)==1, Ky=Ky*ones(size(xm)); end
if size( S,2)==1, S =S *ones(size(xm)); end

Nodes = reshape(1:(Nx+1)*(Ny+1),Ny+1,Nx+1);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
if ~exist('mode','var') || isempty(mode)
    if length(t)>2, fprintf('%s in flat mode.\n',mfilename); end
    Cx     =    dy*(1./dx).*Kx;
    Cyinner=(1./dy)*   dx/2.*Ky;
    Cyouter=(1./dy)*   dx/2.*Ky;
    Csinner=    dy.*   dx/2.*S;
    Csouter=    dy.*   dx/2.*S;
else    
    if length(t)>2, fprintf('%s in radial mode.\n',mfilename); end
    r2=xGr(2:end);     r1=xGr(1:end-1);  rm=0.5*(r2+r1);

    Cx=2*pi*dy*(1./log(r2./r1)).*Kx;
    switch mode
        case 1, % phi is linear with r
            fprintf('linear head assumed with r\n');
            Cyinner=pi*(1./dy)*((2/3*r1.^3+1/3*r2.^3-r1.^2.*r2)./(r2-r1)).*Ky;
            Cyouter=pi*(1./dy)*((2/3*r2.^3+1/3*r1.^3-r2.^2.*r1)./(r2-r1)).*Ky;
            Csinner=pi*    dy *((2/3*r1.^3+1/3*r2.^3-r1.^2.*r2)./(r2-r1)).*S;
            Csouter=pi*    dy *((2/3*r2.^3+1/3*r1.^3-r2.^2.*r1)./(r2-r1)).*S;
        case 2, % phi is logarithmic with r
            fprintf('logarithmic head assumed with r\n');
            LOG=(r2.^2-r1.^2)./log(r2./r1)/2;
            Cyinner=pi*(1./dy)*(LOG - r1.^2).*Ky;
            Cyouter=pi*(1./dy)*(r2.^2 - LOG).*Ky;
            Csinner=pi*    dy *(LOG - r1.^2).*S;
            Csouter=pi*    dy *(r2.^2 - LOG).*S;
        otherwise, % old default, simple
            fprintf('simple vertical conductance\n');
            Cyinner=pi*(1./dy)*(rm.^2-r1.^2).*Ky;
            Cyouter=pi*(1./dy)*(r2.^2-rm.^2).*Ky;
            Csinner=pi*    dy *(rm.^2-r1.^2).*S;
            Csouter=pi*    dy *(r2.^2-rm.^2).*S;
    end
end
Cx=([zeros(size(dx)) ;Cx     ]+[Cx;     zeros(size(dx)) ])/2;
Cy=([zeros(size(dy))  Cyinner]+[Cyouter zeros(size(dy)) ]);
Cs=([zeros(size(dy))  Csinner]+[Csouter zeros(size(dy)) ]);
Cs=([zeros(size(xGr));Cs     ]+[Cs;     zeros(size(xGr))])/2;

Cs=Cs(:);  % storage conductacne when devided by dt*theta

A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:)],...
         (Ny+1)*(Nx+1),(Ny+1)*(Nx+1),5*(Ny+1)*(Nx+1));                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct =Nodes(IBOUND~=0); % active cells including fixed heads
I    =Nodes(IBOUND> 0); % active cells but not fixed heads = cells with heads to be computed
Ifh  =Nodes(IBOUND< 0); % active cells with fixed heads

Phi=NaN((Ny+1)*(Nx+1),Nt+1);  % allocate space to store the entire head matrix
Qt =NaN((Ny+1)*(Nx+1),Nt);    % allocate memory for Qt
Qs =NaN((Ny+1)*(Nx+1),Nt);    % allocate memory for Qs
Qx =NaN((Ny+1),Nx  ,Nt);  % allocate memory for Qx
Qy =NaN( Ny   ,Nx+1,Nt);  % allocate memory for Qy

Phi(IAct,1)=IH(IAct);  % store initial head at Phi(:,:,1)
FQ=FQ(:); 
IH =IH(:);
Fi =IH;            % head computed in this time step
for it=1:Nt
    if isempty(Ifh)
        Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I)                 +Cs(I).*Phi(I,it)/dt(it)/theta); % solve
    else
        Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I)-A(I,Ifh)*IH(Ifh)+Cs(I).*Phi(I,it)/dt(it)/theta); % solve
    end
    Phi(IAct ,it+1)=Fi(IAct)/theta-(1-theta)/theta*Phi(IAct,it);
    Qt (IAct ,it  )=spdiags(Adiag(IAct),0,A(IAct,IAct))*Fi(IAct);
    Qs (IAct ,it  )=-Cs(IAct).*(Phi(IAct,it+1)-Phi(IAct,it))/dt(it);   % Storage in time step m3 for cell
    Qx (:,: ,it   )=-Cx.*diff(reshape(Fi,Ny+1,Nx+1),1,2);	% Flow across horizontal cell faces m3/d per cell
    Qy (:,: ,it   )=-Cy.*diff(reshape(Fi,Ny+1,Nx+1),1,1);   % Flow across vertical cell faces, m3/d per cell
end
Phi=reshape(Phi,(Ny+1),(Nx+1),Nt+1);                   % NaN if inactive
Qt =reshape(Qt ,(Ny+1),(Nx+1),Nt  ); Qt(isnan(Qt))=0;  % 0 if inactive
Qs =reshape(Qs ,(Ny+1),(Nx+1),Nt  ); Qs(isnan(Qs))=0;  % 0 if inactive

Qx(isnan(Qx))=0;
Qy(isnan(Qy))=0;              % 0 if inactive
