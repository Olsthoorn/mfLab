function [Phi,t,Qt,Qx,Qy,Qs]=fdm2tsp(x,y,Period,kx,ky,S,IH,FH,FQ,varargin)
%FDM2TSP a 2D block-centred transient finite difference model with sress periods
%
% Example:
%    [Phi,t,Q,Qx,Qy,Qs]=fdm2tsp(x,y,Period,kx,ky|[],S,IH,FH|[],FQ|[],[,radial])
%
% INPUT:
%  x(Nx+1)         x-coordinate of mesh/grid
%  y(Ny+1)         y-coordinate of mesh/grid
%  Period          Stress periods as a [Nper,3] matrix holding
%  [PerLen, nSTep, Mult]  period length, number of timesteps and time step multiplyer
%  kx(Ny,Nx)       conductivity in x-direction
%  ky(Ny,Nx)       same in y direction, ky=kx if ky==[]
%  S(Ny,Nx)        specifice storage (Ss+Sy)
%  IH(Ny,Nx)       initial head
%  FH --- fixed heads  no fixed heads if []
%     FH can be given as multi layers [Ny,Nx,Nper], one layer per stres period
%        or FH is given in Modflow style
%        one record with th e number of cells to be defined for the next
%        stress period holding [iLay, iRow, iCol, FH]
%
% FQ ---- fixed nodal flows no fixed flows (all zeros) if []
%    FQ can be given as multi-layers [Ny,Nx,Nper], one layer per stress period
%       or fQ is given in Modflow style
%       one record with the number of cells to be defined for the next period
%       holdin g [iLay iRow iCol FH]
% Radial          Arbitrary input caused model to assume axial flow with
%                 x=r and FQ are ring flows while als flows have dimension
%                 L3/T instead of L2/t. For clearness use 'R' or 'Radial'
%                 as input at the function call.
% OUTPUT:
%  Phi(Ny,Nx,Nt+1) computed heads with Phi(Ny,Nx,1) initial heads for t=0
%  t = times of head output
%, Qt(Ny,Nx,Nt)    computed total cell balance during timestep it
%  Qx(Ny,Nx-1,Nt)  hor.  cell face flow in x-direction positive along increasing col number
%  Qy(Ny-1,Nx,Nt)  vert. cell face flow in y-directin  positive along increasing row number
%  Qs(Ny,Nx,Nt)    storage change of node during timestep it
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
% TO 080226 (added inactive cells and true fixed heads, t(1) is now always made zero

% Copyright 2009 Theo Olsthoorn, without any warranty
% under free software foundnation GNU version 3 or later licence
% TO 090328  adding stres periods

% Copyright 2009 Theo Olsthoorn, without any warranty
% under free software foundnation GNU version 3 or later licence

%% implicitness
theta=1;

%% ky
if isempty(ky), ky=kx; end

%% Check stress period input for time
if size(Period,2)~=3
    error('Stress period format is [tend Nstep StpMult] one line for each stress period');
else
    if any(Period(:,1)<=0), error('Stress periods must all be > 0'); end
    if any(Period(:,2)< 1), error('Time steps per stress period must all > 0'); end
    if any(Period(:,3)<=0), error('Time step multiplyers must all be > 0'); end
    Nper=size(Period,1);
end

%% housekeeping grid
x=x(:)'; dx=    diff(x);  Nx=length(dx); xm=0.5*(x(1:end-1)+x(2:end));
y=y(:);  dy=abs(diff(y)); Ny=length(dy);

Nodes = reshape(1:Nx*Ny,Ny,Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

%% Active cells
IAct =Nodes( kx>0 | ky>0 | S>0);               % active cells

%% resistances and conducctances
if nargin>9
    fprintf('%s in radial mode.\n',mfilename);
    RX=(1./dy)*log(x(2:end-1)./xm(1:end-1))./(2*pi*kx(:,1:end-1))+...
       (1./dy)*log(xm(2:end)./x(2:end-1)) ./(2*pi*kx(:,2:end));
    RY=0.5/pi*dy*(1./(x(2:end).^2-x(1:end-1).^2))./ky;
    Cx=1./RX; Cy=1./(RY(1:end-1,:)+RY(2:end,:));
    Cs=pi*dy*(x(2:end).^2-x(1:end-1).^2).*S;
else
    fprintf('%s in flat mode.\n',mfilename);
    RX=0.5*(1./dy)*dx./kx; Cx=1./(RX(:,1:end-1)+RX(:,2:end));
    RY=0.5*dy*(1./dx)./ky; Cy=1./(RY(1:end-1,:)+RY(2:end,:));
    Cs=dy*dx.*S;
end
Cs=Cs(:);  % storage conductacne when devided by dt*theta

%% System matrix
A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         -[Cx(:);Cx(:);Cy(:);Cy(:)],...
         Ny*Nx,Ny*Nx,5*Ny*Nx);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

%% FH, dealing with it depends on the way it was passed into the model Modflow or layer wise
if isempty(FH), FH=NaN(Ny,Nx,Nper); end
if size(FH,2)~=5 && size(FH,1)==Ny && size(FH,2)==Nx && size(FH,3)==Nper
    FH=reshape(FH,Ny*Nx,Nper);
else
    fh=FH;  FH=NaN(Ny*Nx,Nper);
    if any(fh(:,2)~=1), error('Layer Nr in FH must be 1 in this model, [iL iR iC]'); end 
    for iPer=1:Nper
        I=find(fh(:,1)==iPer);
        for i=1:length(I)
            FH(fh(I(i),3),fh(I(i),4),iPer)=fh(I(i),5);
        end
    end
    FH=reshape(FH,Ny*Nx,Nper);
end

%% FQ, dealing with it depends on the way it was passed into the model Modflow or layer wise
if isempty(FQ), FQ=zeros(Ny,Nx,Nper); end
if size(FQ,2)~=5 && size(FQ,1)==Ny && size(FQ,2)==Nx && size(FQ,3)==Nper
    FQ=reshape(FQ,Ny*Nx,Nper);
else
    fq=FQ;  FQ=zeros(Ny,Nx,Nper);
    if any(fq(:,2)~=1), error('Layer Nr in FQ must be 1 in this model, [iL iR iC]'); end 
    for iPer=1:Nper
        I=find(fq(:,1)==iPer);
        for i=1:length(I)
            FQ(fq(I(i),3),fq(I(i),4),iPer)=fq(I(i),5);  % iL=fq(I(2)) is nit used
        end
    end
    FQ=reshape(FQ,Ny*Nx,Nper);
end

%% Initialize output
Nt = sum(Period(:,2));
Phi=NaN(Ny*Nx,Nt+1);  % allocate space to store the entire head matrix
Qt =NaN(Ny*Nx,Nt);    % allocate memory for Qt
Qs =NaN(Ny*Nx,Nt);    % allocate memory for Qs
Qx =NaN(Ny,Nx-1,Nt);  % allocate memory for Qx
Qy =NaN(Ny-1,Nx,Nt);  % allocate memory for Qy
Fi =NaN(Ny,Nx);

%% Solve
k=1;
t=zeros(1,sum(Period(:,2))+1);
for iPer=1:Nper;
    dt=ones(1,Period(iPer,2));   for i=2:Period(iPer,2); dt(i)=dt(i-1)*Period(iPer,3); end; dt=dt*Period(iPer,1)/sum(dt);
    I    =Nodes((kx>0 | ky>0 | S>0) &  isnan(reshape(FH(:,iPer),Ny,Nx))); % active cells but not fixed heads = cells with heads to be computed
    Ifh  =Nodes((kx>0 | ky>0 | S>0) & ~isnan(reshape(FH(:,iPer),Ny,Nx))); % active cells with fixed heads
    if iPer==1
        Phi(IAct,k)=IH(IAct);
        Phi(Ifh ,k)=FH(Ifh,iPer);
    end
    for it=1:length(dt)
        t(k+1)=t(k)+dt(it);
        Fi (I)        =spdiags(Adiag(I)+Cs(I)/dt(it)/theta,0,A(I,I))\(FQ(I,iPer)-A(I,Ifh)*FH(Ifh,iPer)+Cs(I).*Phi(I,k)/dt(it)/theta); % solve
        Phi(IAct ,k+1)=Fi(IAct)/theta-(1-theta)/theta*Phi(IAct,k);
        Qt (IAct ,k  )=spdiags(Adiag(IAct),0,A(IAct,IAct))*Fi(IAct);
        Qs (IAct ,k  )=-Cs(IAct).*(Phi(IAct,k+1)-Phi(IAct,k))/dt(it);   % Storage in time step m3 for cell
        Qx (:,:  ,k  )=-Cx.*diff(reshape(Fi,Ny,Nx),1,2);	% Flow across horizontal cell faces m3/d per cell
        Qy (:,:  ,k  )=-Cy.*diff(reshape(Fi,Ny,Nx),1,1);   % Flow across vertical cell faces, m3/d per cell
        k=k+1;
    end
end
Phi=reshape(Phi,Ny,Nx,Nt+1);                   % NaN if inactive
Qt =reshape(Qt ,Ny,Nx,Nt  ); Qt(isnan(Qt))=0;  % 0 if inactive
Qs =reshape(Qs ,Ny,Nx,Nt  ); Qs(isnan(Qs))=0;  % 0 if inactive
Qx(isnan(Qx))=0; Qy(isnan(Qy))=0;              % 0 if inactive
