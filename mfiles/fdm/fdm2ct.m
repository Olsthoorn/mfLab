function [Phi,Qt,Qx,Qy,Qs]=fdm2ct(xGr,yGr,t,Kx,c,Ky,S,IBOUND,IH,FQ,varargin)
%FDM2CT a 2D block-centred transient finite difference model with semi-confined layers
%
% Example
%    [Phi,Q,Qx,Qy,Qs]=fdm2ct(xGr,yGr,t,Kx,c,Ky,S,IH,IH,FQ [,radial])
%
% INPUTS:
%  xGr(Nx+1)       xGr-coordinate of mesh/grid
%  yGr(Ny+1)       yGr-coordinate of mesh/grid
%  t(Nt+1)         time for output. 0 will be added if necessary
%  Kx(Ny,Nx)       conductivity in x-direction
%  C(Ny-1,Nx)      confining bed resistance
%  Ky(Ny,Nx)       same in y direction, Ky=Kx if Ky==[]
%  S(Ny,Nx)        primary storage (S+Sy)
%  IH(Ny,Nx)       initial head
%  IH(Ny,Nx)       fixed heads (NaN for ordinary points), Q=fixed nodal flows
%  Radial          Arbitrary input caused model to assume axial flow with
%                 xGr=r and FQ are ring flows while als flows have dimension
%                 L3/T instead of L2/t. For clearness use 'R' or 'Radial'
%                 as input at the function call.
% OUTPUTS:
%  Phi(Ny,Nx,Nt+1) computed heads with Phi(Ny,Nx,1) initial heads for t=0
%, Qt(Ny,Nx,Nt)    computed total cell balance during timestep it
%  Qx(Ny,Nx-1,Nt)  hor.  cell face flow in x-direction positive along increasing col number
%  Qy(Ny-1,Nx,Nt)  vert. cell face flow in y-direction positive along increasing row number
%  Qs(Ny,Nx,Nt)    storage change of node during timestep it
%
% See also: fdm2 fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
% TO 080226 (added inactive cells and true fixed heads, t(1) is now always made zero

% Copyright 2010 Theo Olsthoorn, without any warranty
% under free software foundnation GNU version 3 or later licence

theta= 0.67; % implicitness

[xGr,~,xm,~,dx,dy,Nx,Ny]=modelsize(xGr,yGr);

t=unique(t(:));  dt=    diff(t);  Nt=length(dt);

if all(IBOUND)==0, error('IBOUND all zeros (inactive), nothing to compute!'); end

Kx = augment(Kx,Ny,Nx);
Ky = augment(Ky,Ny,Nx);
S  = augment(S ,Ny,Nx);
c  = augment(c ,Ny,Nx);
IH = augment(IH,Ny,Nx);
FQ = augment(FQ,Ny,Nx);

Nodes = reshape(1:Nx*Ny,Ny,Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
if isempty(varargin)
    if length(t)>2, fprintf('%s in flat mode.\n',mfilename); end
    RX=0.5*(1./dy)*    dx ./Kx;
    RY=0.5*    dy *(1./dx)./Ky;
    
    if ~exist('c','var') || isempty(c),
        Rc=zeros(Ny-1,Nx);
    else
        Rc=     c./dx;
    end % confining bed resistance
    
    Cx=1./(RX(:,1:end-1)     +RX(:,2:end));
    
    if Ny==1,
        Cy=[];
    else
        Cy=1./(RY(1:end-1,:)+ Rc +RY(2:end,:));
    end
    
    Cs=(dy*dx).*S;
else    
    if length(t)>2, fprintf('%s in radial mode.\n',mfilename); end
    Area=pi*(xGr(:,2:end).^2-xGr(:,1:end-1).^2);
%    Area(  1)=xGr( 2)^2-xm(     1)^2;  % first cell only part area
%    Area(end)=xm(end)^2-xGr(end-1)^2;  % last  cell only part area

    RX=dy*log(xGr(:,2:end-1)./xm( :,1:end-1))./(2*pi*Kx(:,1:end-1))+...
       dy*log(xm( :,2:end  )./xGr(:,2:end-1))./(2*pi*Kx(:,2:end  ));
    RY=0.5*dy*(1./Area)./Ky;
    
    if ~exist('c','var') || isempty(c),
        Rc=zeros(Ny-1,Nx);
    else
        Rc=c./Area;
    end % confining bed resistance
    
    Cx=1./ RX;
    if Ny==1,
        Cy=[];
    else
        Cy=1./(RY(1:end-1,:)+ Rc +RY(2:end,:));
    end
    Cs=dy*Area.*S;
end
Cs=Cs(:);  % storage conductacne when devided by dt*theta

A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:)],...
         Ny*Nx,Ny*Nx,5*Ny*Nx);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct =Nodes(IBOUND~=0); % active cells
I    =Nodes(IBOUND >0); % active cells but not fixed heads = cells with heads to be computed
Ifh  =Nodes(IBOUND <0); % active cells with fixed heads

Phi=NaN(Ny*Nx,Nt+1);  % allocate space to store the entire head matrix
Qt =NaN(Ny*Nx,Nt);    % allocate memory for Qt
Qs =NaN(Ny*Nx,Nt);    % allocate memory for Qs
Qx =NaN(Ny,Nx-1,Nt);  % allocate memory for Qx
Qy =NaN(Ny-1,Nx,Nt);  % allocate memory for Qy

Phi(IAct,1)=IH(IAct);  % store initial head at Phi(:,:,1)
FQ=FQ(:); 
Fi =IH(:);            % head computed in this time step
if any(I(:))
    for it=1:Nt
        if isempty(Ifh)
            Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I)                 +Cs(I).*Phi(I,it)/dt(it)/theta); % solve
        else
            Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I)-A(I,Ifh)*IH(Ifh)+Cs(I).*Phi(I,it)/dt(it)/theta); % solve
        end
        Phi(IAct ,it+1)=Fi(IAct)/theta-(1-theta)/theta*Phi(IAct,it);
        Qt (IAct ,it  )=spdiags(Adiag(IAct),0,A(IAct,IAct))*Fi(IAct);
        Qs (IAct ,it  )=-Cs(IAct).*(Phi(IAct,it+1)-Phi(IAct,it))/dt(it);   % Storage in time step m3 for cell
        if ~isempty(Cx)
            Qx (:,: ,it   )=-Cx.*diff(reshape(Fi,Ny,Nx),1,2);	% Flow across horizontal cell faces m3/d per cell
        end
        if ~isempty(Qy)
            Qy (:,: ,it   )=-Cy.*diff(reshape(Fi,Ny,Nx),1,1);   % Flow across vertical cell faces, m3/d per cell
        end
    end
end
Phi=reshape(Phi,Ny,Nx,Nt+1);                   % NaN if inactive
Qt =reshape(Qt ,Ny,Nx,Nt  ); Qt(isnan(Qt))=0;  % 0 if inactive
Qs =reshape(Qs ,Ny,Nx,Nt  ); Qs(isnan(Qs))=0;  % 0 if inactive

if ~isempty(Cx), Qx(isnan(Qx))=0; end
if ~isempty(Cy), Qy(isnan(Qy))=0; end              % 0 if inactive

end

function V = augment(V,Ny,Nx)
    if size(V,2)==1
        if size(V,1)==1
            V = V*ones(Ny,Nx);
        else
            V = bsxfun(@times,V,ones(1,Nx));
        end
    end
end

