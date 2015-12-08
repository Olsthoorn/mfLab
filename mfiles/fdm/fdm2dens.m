function [Phi,Q,Psi,Qx,Qy]=fdm2dens(x,y,kx,ky,FH,FQ,gamma)
%FDM2DENS a 2D block-centred steady-state finite difference model with density flow
%
% Example:
%    [Phi,Q,Psi,Qx,Qy]=fdm2d(x,y,kx,ky|[],FH,FQ,gamma [,radial)
%
% x(Nx+1),y(Ny+1) mesh coordinates,kx,ky conductivities
% Nx,NY being the number of cells (not used in the function call)
%
% All heads are point water heads, heads as measured with the local density
% Phif=y+(Phi-y)*rho/rhof=y+(Phi-y)*(1+gamma);
%
% kx(Ny,Nx) and ky(Ny,Nx) cell conductivities, if ky==[] then kx=ky
% FH(Ny,Nx)=fixed heads (NaN for ordinary points)
% FQ(Ny,Nx)=fixed nodal flows
% gamma(Ny,Nx)= relative cell density, rho/rho0
% Radial is an arbirary string: e.g. use 'R' or 'radial'. The presence of a
% string causes fdm2dens to compute axial symmetric flow, interpreting x as r
% and also interpreting Q, Psi and FQ as ring values with dimension [L3/T] instead of [L2/T]
%
% OUTPUTS:
%    Phi,Q computed heads and cell balances (Ny,Nx)
%    Psi is the stream function (Ny,Nx-1)
%    Qx is flowrightface (Ny,Nx-1) postivie in positive x-diretion
%    Qy is flowlowerface (Ny-1,Nx) positive in postivie y-direction
%
% See also: fdm2 fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 090228(add active cells)


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('hellow he\n');

s=sign(y(end)-y(1));
if s<0, y=flipud(y); kx=flipud(kx); ky=flipud(ky); FH=flipud(FH); FQ=flipud(FQ); gamma=flipud(gamma); end

x=x(:)'; Nx=length(x)-1; dx=diff(x); xm=0.5*(x(1:end-1)+x(2:end));
y=y(:);  Ny=length(y)-1; dy=abs(diff(y)); ym=0.5*(y(1:end-1)+y(2:end));

YM=ym*ones(size(xm));  % cell center elevations, all cells
FH=YM+(FH-YM).*(1+gamma);  % convert to fresh-water heads

Nodes = reshape(1:Nx*Ny,Ny,Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

Gam=zeros(Ny,Nx);
if size(gamma,1)~=Ny || size(gamma,2)~=Nx, gamma=Gam; end
Gam(1,:)=dy(1)*(gamma(1,:)+gamma(2,:))/2;
for i=2:Ny
    Gam(i,:)=Gam(i-1,:)+dy(i)*(gamma(i,:)+gamma(i-1,:))/2;
end

% resistances and conducctances
switch nargin
    case 8
        fprintf('Fdm2dens in radial mode.\n');
        RX=(1./dy)*log(x(2:end-1)./xm(1:end-1))./(2*pi*kx(:,1:end-1))+...
           (1./dy)*log(xm(2:end)./x(2:end-1)) ./(2*pi*kx(:,2:end));
        RY=0.5/pi*dy*(1./(x(2:end).^2-x(1:end-1).^2))./ky;
        Cx=1./RX; Cy=1./(RY(1:end-1,:)+RY(2:end,:));
    otherwise
        fprintf('Fdm2dens in flat mode.\n');
        RX=0.5*(1./dy)*dx./kx; Cx=1./(RX(:,1:end-1)+RX(:,2:end));
        RY=0.5*dy*(1./dx)./ky; Cy=1./(RY(1:end-1,:)+RY(2:end,:));
end

% Density gradient effect
b=zeros(size(FQ));
b(1:end-1,:)=b(1:end-1,:)+Cy.*(Gam(2:end,:)-Gam(1:end-1,:));
b(2:end  ,:)=b(2:end  ,:)-Cy.*(Gam(2:end,:)-Gam(1:end-1,:));

%Active and fixed head cells
Iact=find( kx>0 | ky>0);                % the active cells, 2D steadt state
Ifh =find((kx>0 | ky>0) & ~isnan(FH));  % the fixed head cells
I   =find((kx>0 | ky>0) &  isnan(FH));  % active and not fixed head (to be computed)


A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:)],...
         Ny*Nx,Ny*Nx,5*Ny*Nx);                 % System matrix

A=spdiags(-sum(A,2),0,A);                      % add diagonal
% note if diagonal is to be changed, keep it outside !!

Phi=FH;
Q  =zeros(Ny,Nx);

Phi(I)=A(I,I)\(FQ(I)-A(I,Ifh)*FH(Ifh)+b(I)); % solve, including b and d due to density

Q(Iact)=A(Iact,Iact)*Phi(Iact)-b(Iact);
%reshape(spdiags(Adiag,0,A)*Phi(:),Ny,Nx)-b;		% reshape

Phi=YM+(Phi-YM)./(1+gamma);  % convert from fresh water heads to point water heads

Qx=-Cx.*diff(Phi,1,2)                *sign(x(end)-x(1));   % Flow across horizontal cell faces
Qy=-Cy.*(diff(Phi,1,1)+diff(Gam,1,1))*sign(y(end)-y(1));   % Flow across vertical cell faces

if s<0, Phi=flipud(Phi); Q=flipud(Q); Qx=flipud(Qx); Qy=flipud(Qy); end

Psi=[flipud(cumsum(Qx(end:-1:1,:))); zeros(size(Qx(1,:)))];
