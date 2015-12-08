function [xfm,zfm,gamma]=getinterface(x,z,delta,xf,zf)
% GETINTERFACE Yields coordinates of an interface
%
% USAGE:
%   [xf,zf,gamma]=getinterface(x,z,delta,xf,zf)
%
%   x(Nx+1),z(Nz+1) grid coordinates, any value, usually (rho-rho0)/rho
%   xf,zf optional interface coordinates, if omitted the interface must be
%   clicked in with the mouse.
%
%   The area above the interface gets density zero,
%   the area below it gets density delta
%
%   gamma is the matrix of densities at the center of the horizontal cell faces.
%   If the line is more or less a closed polygon, it will closed automatically and the
%   encirculated area gets the density. If density (rho-rho0)/rho < 0  the
%   encircled area is lighter than the surrounding fluid.
%
% TO 070521 070616

% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


s=sign(z(end)-z(1)); if s<0, z=flipud(z); end

Nx=length(x)-1; dx=diff(x); %xm=0.5*(x(1:end-1)+x(2:end));
Nz=length(z)-1; dz=diff(z); %zm=0.5*(z(1:end-1)+z(2:end));

if nargin<4
    xf=([]); zf=([]); 
    i=1;
%    while 1
%        [xf(i),zf(i),button]=ginput(1);
%        if button~=1
%             xf(end)=[]; zf(end)=[];
%            break;
%        elseif i>1
%            hl=line(xf(end-1:end),zf(end-1:end),'color','b');
%        end
%        i=i+1;
%    end
    [xf,zf]=getline(gca); % from the images toolbox
else
   hl=line(xf(end-1:end),zf(end-1:end),'color','b');
end

if strcmp(get(gca,'xscale'),'log')
    xf=log10(xf);
    x =log10(x);
    dx=diff(x);
end

% remove same coordinates
d=abs(diff(xf))+abs(diff(zf)); I=find(d<=eps);
for i=length(I):-1:1, xf(I(i))=[]; zf(I(i))=[]; end


r=sqrt((xf-xf(1)).^2+(zf-zf(1)).^2);
if r(end)<0.1*max(r)
    xf(end)=xf(1); zf(end)=zf(1);
    if strcmp(get(gca,'xscale'),'log')
        set(hl,'xdata',10.^xf(end-1:end),'ydata',zf(end-1:end));
    else
        set(hl,'xdata',    xf(end-1:end),'ydata',zf(end-1:end));
    end
    drawnow;
end

xfm=([]); zfm=([]); ifm=1;  % detailed interface points

A=zeros(Nz,Nx);     % computed cell area below interface
Ao=dz*dx;           % total cell area

% ==== intersect the interface with the grid

Np=length(xf);
if Np>1
    for ip=1:Np-1  % all vectors in interface points
        x1=xf(ip); x2=xf(ip+1); dx=x2-x1;
        z1=zf(ip); z2=zf(ip+1); dz=z2-z1;

        if abs(dx)>eps        % vector is nor vertical
            lambda=(x([1,end])-x1)/dx;
            % clip on left and right model mesh borders
            if lambda(2)<1+eps && lambda(2)>-eps && lambda(1)<lambda(2) 
                x2=x1+lambda(2)*dx;
                z2=z1+lambda(2)*dz;
            end
            if lambda(1)<1+eps && lambda(1)>-eps && lambda(2)>lambda(1)
                x1=x1+lambda(1)*dx;
                z1=z1+lambda(1)*dz;
            end
            if lambda(1)<1+eps && lambda(1)>-eps && lambda(2)<lambda(1) 
                x2=x1+lambda(1)*dx;
                z2=z1+lambda(1)*dz;
            end
            if lambda(2)<1+eps && lambda(2)>-eps && lambda(1)>lambda(2)
                x1=x1+lambda(2)*dx;
                z1=z1+lambda(2)*dz;
            end
            dx=x2-x1;
            dz=z2-z1;
        end
        
        if abs(dx)>eps
            lambda=(x([1,end])-x1)/dx;
        elseif  x1>=x(1) && x1<=x(end)   % interface vector is vertical but inside left and right model boundary
            lambda=[0 1];              % dummy
        else
            lambda=[NaN NaN];          % vertical interface outside model, ignore
        end
        
        if (lambda(1)<=eps && lambda(2)>=1-eps) || (lambda(1)>=1-eps && lambda(2)<=eps) % vector is within odel boundary
            lambda=[];
            if abs(dx)>eps, lambda=[lambda (x-x1) /dx]; end
            if abs(dz)>eps, lambda=[lambda (z-z1)'/dz]; end
            lambda=sort(lambda(lambda>eps & lambda<1-eps));
            if ~isempty(lambda)                       % if none then what, use the begin and end point directly
                    xp=[x1, x1+lambda*dx, x2];
                    zp=[z1, z1+lambda*dz, z2];
            else
                xp=[x1 x2];
                zp=[z1 z2];
            end
            ddx=diff(xp);
            for ipp=1:length(ddx);
                xmm=0.5*(xp(ipp+1)+xp(ipp));
                zmm=0.5*(zp(ipp+1)+zp(ipp)); 

                xfm(ifm)=xmm;
                zfm(ifm)=zmm;
                ifm=ifm+1;

                ic=find(x(1:Nx)<xmm,1,'last');
                J=find(z(1:Nz)<zmm);
                for jj=1:length(J)
                    jc=J(jj);
                    A(jc,ic)=A(jc,ic)+ddx(ipp)*min(zmm-z(jc),z(jc+1)-z(jc));
                end
            end
        end
    end
end

gamma=delta*A./Ao;
if s<0; gamma=flipud(gamma); end

if strcmp(get(gca,'xscale'),'log')
    xfm=10.^xfm;
end
