function [Phi,Q,Psi,Qx,Qy]=fdm2NC(xGr,yGr,Kx,Ky,IBOUND,IH,FQ,varargin)
%FDM2NC a 2D node-centred steady-state finite difference model
%
% Example:
%    [Phi,Q,Psi,Qx,Qy]=fdm2d(xGr,yGr,Kx,Ky,IBOUND,IH,FQ [,radial])
%
% INPUTS:
% xGr,yGr mesh coordinates, Kx,Ky conductivities, if Ky==[], Ky=Kx
% xGr must be a hor  vector or a grid
% yGr must be a vert vector or a grid
% Kx=horizontal conductivities
% Kv=vertical conductivities
% IBOUND array showing which nodes to compute (>0), which are fixed (<0)
% and which are inactive (0)
% IH= initial heads for nodes
%
% OUTPUTS:
% FQ= node injections
% Phi,Q computed heads and cell balances
% Qx(Ny-1,Nx-2) horizontal cell face flow positive in positive xGr direction
% Qy(Ny-2,Nx-1) vertial    cell face flow, postive in positive yGr direction
% Psi(Ny,Nx-2)  stream function
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 090314 101130 120111

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if ~isempty(varargin)
    mode=varargin{1};
else
    mode=2;
end
    
[xGr,yGr,~,~,dx,dy,Nx,Ny]=modelsize(xGr,yGr);

if numel(Kx)==1, Kx=ones(Ny,Nx)*Kx; elseif numel(Kx(1,:))==1, Kx=Kx*ones(1,Nx); end
if numel(Ky)==1, Ky=ones(Ny,Nx)*Ky; elseif numel(Ky(1,:))==1, Ky=Ky*ones(1,Nx); end

NOD=(Nx+1)*(Ny+1);

Nodes = reshape(1:NOD,Ny+1,Nx+1);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
if isempty(varargin)
    fprintf('%s in flat mode.\n',mfilename)
    Cx     =Kx.*dy*(1./dx);
    Cyinner=Ky.*(1./dy)*dx/2;
    Cyouter=Ky.*(1./dy)*dx/2;
else
    fprintf('%s in axisymmetric mode.\n',mfilename)
    r2=xGr(2:end);
    r1=xGr(1:end-1);
    Cx=2*pi*dy*(1./log(r2./r1)).*Kx;
    switch mode
        case 1, % phi is linear with r
            fprintf('linear head assumed with r\n');
            Cyinner=pi*(1./dy)*((2/3*r1.^3+1/3*r2.^3-r1.^2.*r2)./(r2-r1)).*Ky;
            Cyouter=pi*(1./dy)*((2/3*r2.^3+1/3*r1.^3-r2.^2.*r1)./(r2-r1)).*Ky;
        case 2, % phi is logarithmic with r
            fprintf('logarithmic head assumed with r\n');
            LOG=(r2.^2-r1.^2)./log(r2./r1)/2;
            Cyinner=pi*(1./dy)*(LOG - r1.^2).*Ky;
            Cyouter=pi*(1./dy)*(r2.^2 - LOG).*Ky;
        otherwise, % old default, simple
            fprintf('simple vertical conductance\n');
            rm=0.5*(r2+r1);
            Cyinner=pi*(1./dy)*(rm.^2-r1.^2).*Ky;
            Cyouter=pi*(1./dy)*(r2.^2-rm.^2).*Ky;
    end
end
Cx=0.5*([zeros(size(dx));Cx]+...
        [Cx;zeros(size(dx))]);
Cy=[zeros(size(dy)) Cyinner]+...
   [Cyouter zeros(size(dy))];

A=sparse([IE(:); IW(:); IN(:); IS(:)],...
         [IW(:); IE(:); IS(:); IN(:)],...
        -[Cx(:); Cx(:); Cy(:); Cy(:)],...
         (Ny+1)*(Nx+1),(Ny+1)*(Nx+1),5*(Ny+1)*(Nx+1));   % System matrix
Adiag= -sum(A,2);                % Main diagonal

IAct=Nodes(IBOUND~=0);

I   = Nodes(IBOUND >0);
Ifh = Nodes(IBOUND <0);
IAct= Nodes(IBOUND~=0);

Phi=IH;
Q  =zeros(size(IH));

if size(IH,1)==1, % i.e. a one layer horizontal model    
    Phi(I)=spdiags(Adiag(I),0,A(I,I))\(FQ(I)'-A(I,Ifh)*IH(Ifh)'); % solve
    Q(IAct)=spdiags(Adiag(IAct),0,A(IAct,IAct))*Phi(IAct)';		% reshape
else
    Phi(I)=spdiags(Adiag(I),0,A(I,I))\(FQ(I)-A(I,Ifh)*IH(Ifh)); % solve
    Q(IAct)=spdiags(Adiag(IAct),0,A(IAct,IAct))*Phi(IAct);		% reshape
end
    
if size(IH',1)==1, Phi=Phi'; end

Qx=-Cx.*diff(Phi,1,2)*sign(xGr(1,end)-xGr(1,1));   % Flow across vertical   cell faces
Qy=-Cy.*diff(Phi,1,1)*sign(yGr(end,1)-yGr(1,1));   % Flow across horizontal cell faces

Qx(isnan(Qx))=0;
Qy(isnan(Qy))=0;

Psi= [flipud(cumsum(flipud(Qx)));zeros(size(Qx(1,:)))];
