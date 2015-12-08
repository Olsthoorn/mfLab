function [Phi,Qt,Qx,Qy,Qs]=fdm2ct(gr,t,Kx,Ky,Ss,IBOUND,IH,FQ,varargin)
%FDM2CT a 2D block-centred transient finite difference model
%
% Example
%    [Phi,Q,Qx,Qy,Qs]=fdm2ct(gr,t,Kx,Ky,Ss,IH,IH,FQ [,radial])
%
% INPUTS:
%  gr              grid2DObj
%  t(Nt+1)         time for output. 0 will be added if necessary
%  Kx(Ny,Nx)       conductivity in x-direction
%  Ky(Ny,Nx)       same in y direction, Ky=Kx if Ky==[]
%  Ss(Ny,Nx)       Specific storage coefficient
%  IH(Ny,Nx)       initial head
%  IH(Ny,Nx)       fixed heads (NaN for ordinary points), Q=fixed nodal flows
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

t=unique(t(:));  dt=    diff(t);  Nt=length(dt);

if all(IBOUND)==0, error('IBOUND all zeros (inactive), nothing to compute!'); end

Kx = augment(Kx,gr.Ny,gr.Nx);
Ky = augment(Ky,gr.Ny,gr.Nx);
Ss = augment(Ss ,gr.Ny,gr.Nx);
IH = augment(IH,gr.Ny,gr.Nx);
FQ = augment(FQ,gr.Ny,gr.Nx);

Nodes = reshape(1:gr.Nod,gr.Ny,gr.Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
if ~gr.AXIAL
    if ~gr.AXIAL, fprintf('%s in flat mode.\n',mfilename); end
    RX=0.5*(1./gr.dy)*    gr.dx ./Kx./gr.dZ;
    RY=0.5*    gr.dy *(1./gr.dx)./Ky./gr.dZ;
    
    Cx=1./(RX(:,1:end-1)     +RX(:,2:end));
    
    if Ny==1,
        Cy=[];
    else
        Cy=1./(RY(1:end-1,:)+ Rc +RY(2:end,:));
    end
    
    Cs=gr.Vol.*Ss/theta;
else    
    if gr.AXIAL, fprintf('%s in radial mode.\n',mfilename); end

    RX= bsxfun(@rdivide, log(gr.xGr(:,2:end-1)./gr.xm( :,1:end-1)),2*pi*Kx(:,1:end-1).* gr.dY(:,2:end))+...
        bsxfun(@rdivide, log(gr.xm( :,2:end  )./gr.xGr(:,2:end-1)),2*pi*Kx(:,2:end  ).* gr.dY(:,2:end));
    RY= bsxfun(@rdivide,0.5*gr.dY./Ky,gr.Area);
    
    Cx=1./ RX;
    if gr.Ny==1,
        Cy=[];
    else
        Cy=1./(RY(1:end-1,:)+ RY(2:end,:));
    end
    Cs=gr.Vol.*Ss/theta;
end
Cs=Cs(:);  % storage conductacne when devided by dt*theta

A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:)],...
         gr.Nod,gr.Nod,5*gr.Nod);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct =Nodes(IBOUND~=0); % active cells
I    =Nodes(IBOUND >0); % active cells but not fixed heads = cells with heads to be computed
Ifh  =Nodes(IBOUND <0); % active cells with fixed heads

Phi=NaN(gr.Nod,Nt+1);  % allocate space to store the entire head matrix
Qt =NaN(gr.Nod,Nt);    % allocate memory for Qt
Qs =NaN(gr.Nod,Nt);    % allocate memory for Qs
Qx =NaN(gr.Ny,gr.Nx-1,Nt);  % allocate memory for Qx
Qy =NaN(gr.Ny-1,gr.Nx,Nt);  % allocate memory for Qy

Phi(IAct,1)=IH(IAct);  % store initial head at Phi(:,:,1)
FQ=FQ(:); 
Fi =IH(:);            % head computed in this time step
if any(I(:))
    for it=1:Nt
        if isempty(Ifh)
            Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it),0,A(I,I))\(FQ(I)                 +Cs(I).*Phi(I,it)/dt(it)); % solve
        else
            Fi (I) = spdiags(Adiag(I)+Cs(I)/dt(it),0,A(I,I))\(FQ(I)-A(I,Ifh)*IH(Ifh)+Cs(I).*Phi(I,it)/dt(it)); % solve
        end
        Phi(IAct ,it+1)=Fi(IAct)/theta-(1-theta)/theta*Phi(IAct,it);
        Qt (IAct ,it  )=spdiags(Adiag(IAct),0,A(IAct,IAct))*Fi(IAct);
        Qs (IAct ,it  )=-Cs(IAct).*(Phi(IAct,it+1)-Phi(IAct,it))/dt(it);   % Storage in time step m3 for cell
        if ~isempty(Cx)
            Qx (:,: ,it   )=-Cx.*diff(reshape(Fi,gr.size),1,2);	% Flow across horizontal cell faces m3/d per cell
        end
        if ~isempty(Qy)
            Qy (:,: ,it   )=-Cy.*diff(reshape(Fi,gr.size),1,1);   % Flow across vertical cell faces, m3/d per cell
        end
    end
end
Phi=reshape(Phi,[gr.size, Nt+1]);                   % NaN if inactive
Qt =reshape(Qt ,[gr.size ,Nt  ]); Qt(isnan(Qt))=0;  % 0 if inactive
Qs =reshape(Qs ,[gr.size, Nt  ]); Qs(isnan(Qs))=0;  % 0 if inactive

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

