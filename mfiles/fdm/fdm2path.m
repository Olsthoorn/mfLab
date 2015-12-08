function [XP, YP, TP]=fdm2path(x,y,DZ,Q,Qx,Qy,por,T,markers,XStart,YStart)
%FDM2PATH 2D particle tracking.
%
% Example:
%    [XP YP TP]=fdm2path(x,y,[DZ|[]],Q,Qx,Qy,por,T,markers [,XStart,YStart])
%
% To use this function:
%    generate a 2D steady-state model, launch this function using its produced
%    matrixes Q,Qx,Qy and other necessary parameters DX and por
%    if no starting values are used you must click on an existing picture and
%    the flow path will be immedidately drawn with markers at the given time points in T
%    Repeat this for more lines. Click te right hand button to stop
%    Type fdm2path for selftest and demo
%
% INPUT:
%    x y mesh coordinates (Nx+1) (Ny+1) Not the cell centers.
%    DZ thickness of the cells (if empty or character string radial mode is used).
%    So use 'radial' or so at this position to de particle tracking in axially symmetric
%    with x=r. Then the first r values must be>0.
%    Q Qx Qy output of fdm2 (steady state only)
%    por matrix of porosities (scalar is enough)
%    T time points where markers are desired, A marker at t=0 will always be placed
%    Use negative times vor backward tracing
%    markers is a series of markers for the consecutive times
%    e.g.   '>+o*.xsdph^v<'. The series will be repeated if necessary.
%    XP YP TP coordintates of flow paths, there is a [NaN NaN NaN] between
%    consecutive tracks if muliple starting points or clicks are used.
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 070424 070501


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0; selftest; return; end

if nargin<9
    markers='+o*.xsdph^v<>';
end
Lm=length(markers);

Nx=length(x)-1; Ny=length(y)-1;

if isscalar(DZ),  DZ =DZ *ones(Ny,Nx); end
if isscalar(por), por=por*ones(Ny,Nx); end

%first make sure the positive direction of the grid is aligned with the positive x and y directions
if sign(x(end)-x(1))<0, x=fliplr(x); Q=fliplr(Q); Qx=fliplr(Qx); Qy=fliplr(Qy); DZ=fliplr(DZ); por=fliplr(por); end
if sign(y(end)-y(1))<0, y=flipud(y); Q=flipud(Q); Qx=flipud(Qx); Qy=flipud(Qy); DZ=flipud(DZ); por=flipud(por); end

dx=diff(x); dy=diff(y);

% then check which cell are sinks

if T(end)<T(1)  % if times negative then track particles backward in tme
    Qx=-Qx;
    Qy=-Qy;
    T=-T;
end

sinkfrac=0.25;

Qy(isnan(Qy))=0;
Qx(isnan(Qx))=0;
Q(isnan( Q ))=0;

Qthrough=zeros(size(Q));
if ~isempty(Qx)
    Qthrough=Qthrough+[zeros(size(Qx(:,1))),abs(Qx)]+[abs(Qx),zeros(size(Qx(:,1)))];
end
if ~isempty(Qy)
    Qthrough=Qthrough+[zeros(size(Qy(1,:)));abs(Qy)]+[abs(Qy);zeros(size(Qy(1,:)))];
end
sink= Q < -sinkfrac*Qthrough;
%figure; spy(sink)
if isempty(DZ) || ischar(DZ) % then the flow is axially symmetric
    fprintf('Fdmpath in radial mode.\n')
    if ~isempty(Qx)
        A=dy*2*pi*x;
        vx2=[Qx, zeros(size(Qx(:,1)))]./(A(:,2:end)  .*por);
        vx1=[zeros(size(Qx(:,1))), Qx]./(A(:,1:end-1).*por);
        ax=(vx2-vx1)./(ones(size(Qx(:,1)))*diff(x));
    end
    if ~isempty(Qy)
        A=ones(size(dy))*(pi*(x(2:end).^2-x(1:end-1).^2));
        vy2=[Qy; zeros(size(Qy(1,:)))]./(A.*por);
        vy1=[zeros(size(Qy(1,:))); Qy]./(A.*por);
        ay=(vy2-vy1)./(diff(y)*ones(size(Qy(1,:))));
    end
else
    fprintf('Fdmpath in flat mode.\n')
    if ~isempty(Qx)
        vx2=[Qx, zeros(size(Qx(:,1)))]./((dy*ones(size(dx))).*por.*DZ);
        vx1=[zeros(size(Qx(:,1))), Qx]./((dy*ones(size(dx))).*por.*DZ);
        ax=(vx2-vx1)./(ones(size(Qx(:,1)))*diff(x));
    end
    if ~isempty(Qy)
        vy2=[Qy; zeros(size(Qy(1,:)))]./((ones(size(dy))*dx).*por.*DZ);
        vy1=[zeros(size(Qy(1,:))); Qy]./((ones(size(dy))*dx).*por.*DZ);
        ay=(vy2-vy1)./(diff(y)*ones(size(Qy(1,:))));
    end
end

XP=([]); YP=([]); TP=([]); j=1;

% startpoints must be inside model
%Iout=find( XStart<min(x) | XStart>max(x) | YStart<min(y) | YStart>max(y) );
%XStart(Iout)=[];
%YStart(Iout)=[];

while 1
    if exist('XStart','var') && exist('YStart','var')
        if j>length(XStart), break; end
        Xp=XStart(j); Yp=YStart(j);  % get starting points for stream lines
        j=j+1;
    else
        [Xp, Yp, button]=ginput(1);  if button~=1; break; end % get starting points for stream lines
    end
    DT=diff(T(:)); if T(1)~=0, DT=[T(1);DT]; end

    for ip=1:length(Xp);
        xp=Xp(ip); yp=Yp(ip); t=T(1);
        XP=[XP;NaN;xp];
        YP=[YP;NaN;yp];
        TP=[TP;NaN; t];
        iLast=length(TP); % to later plot only this  line

        ic=find(x<xp,1,'last'); if isempty(ic) || ic==length(x), break; end
        jc=find(y<yp,1,'last'); if isempty(jc) || jc==length(y), break; end

        line(xp,yp,'marker',markers(1)); hold on;  % initial marker
        
        for idt=1:length(DT);
            dt=DT(idt);
            while dt>0
                if isempty(Qx)
                    dic=0; dtx=dt;
                else
                    [xpN,dic,dtx]=postime(xp,x(ic),x(ic+1),vx1(jc,ic),vx2(jc,ic),ax(jc,ic),dt);
                end
                if isempty(Qy)
                    djc=0; dty=dt;
                else
                    [ypN,djc,dty]=postime(yp,y(jc),y(jc+1),vy1(jc,ic),vy2(jc,ic),ay(jc,ic),dt);
                end

                [ddt,i]=min([dtx,dty]);

                switch i
                    case 1
                        if ~isempty(Qy)
                            xp=xpN;
                            yp=pos(yp,y(jc),vy1(jc,ic),ay(jc,ic),ddt);
                        end
                        ic=ic+dic;
                    case 2
                        if ~isempty(Qx)
                            xp=pos(xp,x(ic),vx1(jc,ic),ax(jc,ic),ddt);
                            yp=ypN;
                        end
                        jc=jc+djc;
                end

                dt=dt-ddt; t=t+ddt;
                XP=[XP;xp]; YP=[YP;yp]; TP=[TP;t];
                if length(XP)>20000; break; end
              
              if dt==0
                  m=mod(idt+1,Lm); if m==0, m=Lm; end   % the +1 because the first marker is the initial one
                  line(xp,yp,'marker',markers(m)); hold on;
              end

              if sink(jc,ic);
                  break;  % from while
              end

            end

            if sink(jc,ic);
              break; % from for 
            end
        end
        line(XP(iLast:end),YP(iLast:end),'color','g');
    end
end
XP=[XP;NaN]; YP=[YP;NaN]; TP=[TP;NaN];
end

function [xp,dic,dt]=postime(xp,x1,x2,v1,v2,ax,Dt)
EPS=1e-6;

v=v1+ax*(xp-x1);
if abs(v)<EPS
    % ic=ic
    dt=Dt;  % immediately jumpt to end of time step
    dic=0;
    return; % x remains same location
end

if v<0  % point moves to face at left side
    if abs(ax)<EPS  % v will be constant
        dt=(x1-xp)/v;
        if dt>Dt
            dt=Dt;
            xp=xp+v*dt;
            dic=0;
        else
            xp=x1;
            dic=-1;
        end
    elseif v1>=0           % point will never reach left face
        dt=Dt;         % immediately jump to end of time step
        xp=pos(xp,x1,v1,ax,dt); % compute position at Dt
        dic=0; % ic=ic
    else
        dt=tim(xp,x1,x1,v1,ax);
        if dt>Dt
            dt=Dt;
            xp=pos(xp,x1,v1,ax,dt);
            dic=0;
        else
            xp=x1;
            dic=-1;
        end
    end
end

if v>0
    if abs(ax)<EPS
        dt=(x2-xp)/v;
        if dt>Dt
            dt=Dt;
            dic=0;
            xp=xp+dt*v;
        else
            xp=x2;
            dic=+1;
        end
    elseif v2<=0
        dt=Dt;
        xp=pos(xp,x1,v1,ax,dt);
        dic=0;
    else
        dt=tim(xp,x2,x1,v1,ax);  % CHECK
        if dt>Dt
            dt=Dt;
            xp=pos(xp,x1,v1,ax,dt);
            dic=0;
        else
            xp=x2;
            dic=+1;
        end
    end
end
end

function xp=pos(xstart,x1,v1,ax,dt)
EPS=1e-6;
if abs(ax)<EPS
    vx=v1+ax*(xstart-x1);
    xp=xstart+vx*dt;
else         
    xp=x1+(v1/ax+(xstart-x1))*exp(ax*dt)-v1/ax;
end
end
    
function dt=tim(xstart,xtarget,x1,v1,ax)
    dt=1/ax*log((v1+ax*(xtarget-x1))/(v1+ax*(xstart-x1)));
end

function selftest
    eval('help fdm2path')
    clear all; close all
    y=linspace(-2500,2500,22);
    x=linspace(-2500,2500,22);

    [x,y,xm,ym,dx,dy,Nx,Ny]=modelsize(x,y);

    DZ=50;
    k =10;
    n=0.001;

    kx= ones(Ny,Nx)*k*DZ; ky=kx;
    FH=zeros(Ny,Nx)*NaN; FH(:,[1,end])=0;  FH([1,end],:)=0;
    FQ=n*dy*dx;
    [Phi,Q,Qx,Qy]=fdm2(x,y,kx,ky,FH,FQ);
    contour(xm,ym,Phi); hold on

    %Track particles
    por=0.35; DZ=50;
    t=[60 365 3650 25*365 100*365];
    [XP,YP,TP]=fdm2path(x,y,DZ,Q,Qx,Qy,por,t,'...p...p...p');
end