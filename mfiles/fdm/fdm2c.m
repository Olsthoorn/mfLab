function [Phi,Q,Psi,Qx,Qy]=fdm2c(xGr,yGr,kx,c,ky,IBOUND,IH,FQ,varargin)
%FDM2C a 2D block-centred steady-state finite difference model with semi-confined
% resistance between the model rows (can be seen as layers in XSection)
%
% Example
%    [Phi,Q,Psi,Qx,Qy]=fdm2d(xGr,yGr,kx,[c],ky,IH,FQ [,radial])
%
% Inputs:
%    xGr,yGr mesh coordinates, kx,ky conductivities, if ky==[], ky=kx
%      xGr must be a hor  vector or a grid
%      yGr must be a vert vector or a grid
%    kx=horizontal conductivities
%    c =confining bed resistances on top of the aquifers (c(1) top)
%      Therefore, the size of c(Ny-2,Nx) is is interpreted as such
%    Kv=vertical conductivities
%    IH=fixed heads (NaN for ordinary points), Q=fixed nodal flows
% Outputs:
%    Phi,Q computed heads and cell balances
%    Qx(Ny-1,Nx-2) horizontal cell face flow positive in positive xGr direction
%    Qy(Ny-2,Nx-1) vertial    cell face flow, postive in positive yGr direction
%    Psi(Ny,Nx-2)  stream function
%
% See also: fmd2 fdm2c fdm2ct fdm3 fdm3t
%
%
% TO 991017  TO 000530 001026 070414 090314 101130

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

[xGr,yGr,xm,~,dx,dy,Nx,Ny]=modelsize(xGr,yGr);

if isempty(c), c=zeros(Ny-1,1); end

if size(c,1)~=Ny-1,
    error('Number of resistant layers (%d) must be Ny=1=%d',size(c,2),Ny);
end

Nodes = reshape(1:Nx*Ny,Ny,Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
if isempty(varargin)
    fprintf('%s in flat mode.\n',mfilename)
    RX=0.5*bsxfun(@rdivide,dx,dy)./kx;
    RY=0.5*bsxfun(@rdivide,dy,dx)./ky;

    if ~exist('c','var') || isempty(c),
        Rc = zeros(Ny-1,Nx);
    else
        Rc = bsxfun(@rdivide,c,dx);
    end % confining bed resistance

    Cx=1./(RX(:,1:end-1)   +RX(:,2:end));
    Cy=1./(RY(1:end-1,:)+Rc+RY(2:end,:));
else
    fprintf('%s in axisymmetric mode.\n',mfilename)
    Area=pi*(xGr(:,2:end).^2-xGr(:,1:end-1).^2);

    if size(dy,2)==1
        dy = dy*ones(size(dx));
    end
    
    RX=(bsxfun(@rdivide,log(xGr(:,2:end-1)./xm( :,1:end-1)),(2*pi*kx(:,1:end-1).*dy(:,1:end-1)))+...
        bsxfun(@rdivide,log(xm( :,2:end  )./xGr(:,2:end-1)),(2*pi*kx(:,2:end  ).*dy(:,2:end  )))...
    );

    RY=bsxfun(@rdivide,0.5*dy./ky,Area);

    if ~exist('c','var') || isempty(c),
        Rc = zeros(Ny-1,Nx);
    else
        Rc = bsxfun(@rdivide,c,Area);
    end % confining bed resistance

    Cx=1./ RX;
    Cy=1./(RY(1:end-1,:)+Rc+RY(2:end,:));
end

A=sparse([IE(:); IW(:); IN(:); IS(:)],...
         [IW(:); IE(:); IS(:); IN(:)],...
        -[Cx(:); Cx(:); Cy(:); Cy(:)],...
         Ny*Nx,Ny*Nx,5*Ny*Nx);   % System matrix
Adiag= -sum(A,2);                % Main diagonal

Inact  = IBOUND(:)==0;
Iact   = IBOUND(:) >0;
Ifh    = IBOUND(:) <0;

IH=IH(:); Phi=IH(:); FQ=FQ(:); Q=FQ(:);

Phi(Iact)=spdiags(Adiag(Iact),0,A(Iact,Iact))\(FQ(Iact)-A(Iact,Ifh)*IH(Ifh)); % solve

Q(~Inact)=spdiags(Adiag(~Inact),0,A(~Inact,~Inact))*Phi(~Inact);		% reshape

Phi = reshape(Phi,size(IBOUND));
Q   = reshape(Q  ,size(IBOUND));

Qx=-Cx.*diff(Phi,1,2)*sign(xGr(1,end)-xGr(1,1));   % Flow across vertical   cell faces
Qy=-Cy.*diff(Phi,1,1)*sign(yGr(end,1)-yGr(1,1));   % Flow across horizontal cell faces

Qx(isnan(Qx))=0;
Qy(isnan(Qy))=0;

Psi= [flipud(cumsum(flipud(Qx)));zeros(size(Qx(1,:)))];

end
