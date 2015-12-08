function [A,xfm,yfm,wfm,Rshore,Lshore]=polyXgrid(x,y,xyw,Code)

% POLYXGRID: New to be made function [A,xfm,yfm,wfm,RShore,LShore]=polyXgrid(x,y,[ []|[xyw], [L]])
%
% TOFIX:  still under concstruction   
%
% USAGE:
%   New to be made function [A,xfm,yfm,wfm,RShore,LShore]=polyXgrid(x,y,[ []|[xyw], [L]])
%   where xyw=[xf,yf,wf] is a line with optional given width (wf).
%   where x,y is the grid
%   xfm yfm wfm are the intersections with the grid including the points in xyz
%   A is a matrix of size Ny,Nx:
%   if wf is not given in xyw, then A is the area of the cells above xf yf
%   if line is left to right or negative area if right to left.
%   The interior of a closed contour is the surface of the cells in the encicled aera if counterclockwise
%   of negative of clockwise.
%   If w is given, A is the area between the righ and left shore defined by
%   xf,yf (heart line) and width wf.
%   Rshore are the coordinates of the right shore
%   Lshore are the coordinates of the left  shore
%   If last argument is used, A with hold the length of the line in each cell
%    multiplied by the width in the cells if wf is given else wf is 1, which
%   yields the length. Calling examples are
%
%   polyxgrid;                                    % Selftest
%   [A,xfm,yfm,wfm]=polyxgrid(x,y);               % Areas, digitiz contour on screen
%   [A,xfm,yfm,wfm]=polyxgrid(x,y,[],'L');        % Lengths, digitize line on screen
%   [A,xfm,yfm,wfm]=polyxgrid(x,y,[xf,yf]);       % Areas, for given line
%   [A,xfm,yfm,wfm,Rshore,LShore]=polyxgrid(x,y,[xf,yf,wf]);    % Areas of stream
%   [A,xfm,yfm,wfm]=polyxgrid(x,y,[xf,yf],'L');   % Length of line in cells
%   [A,xfm,yfm,wfm]=polyxgrid(x,y,[xf,yf,wf],'L'); % Length times width of line in cells
%
% TO 070521 070616 090624


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0; selftest; return; end

%% Grid bookkeeping
Nx=length(x)-1;
Ny=length(y)-1; dy=abs(diff(y));

%s=sign(y(end)-y(1)); if s>0, y=flipud(y); end

%% default
mode =1;  % compute cell area above given line

%% Split input
if nargin>2
    if size(xyw,2)<2,
        error('nodes input must be [xf yf] or [xf yf wf]');
    else
        xf=xyw(:,1); yf=xyw(:,2);
        if size(xyw,2)>2,
            wf=xyw(:,3);
        else
            wf=ones(size(xf));
        end
    end
end

%% in case no polyline given, get it from screen
if nargin<3 || isempty(xyw)
    [xf,yf]=getline(gca); % from the image toolbox
    wf=ones(size(xf));
end

hold on
plot(xf,yf,'--bs');

%% Close automatially if start and endpoints are less than 3 pixels apart
Fpos=get(gcf,'position'); Fw=Fpos(3); Fh=Fpos(4);
Apos=get(gca,'position'); Aw=Apos(3); Ah=Apos(4);
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
dx2pix=2/Fw*Aw*(xlim(2)-xlim(1));
dy2pix=2/Fh*Ah*(ylim(2)-ylim(1));

if abs(xf(end)-xf(1))<=dx2pix && abs(yf(end)-yf(1))<dy2pix
    xf(end)=xf(1);
    yf(end)=yf(1);
end

%% Remove points that are double
d=abs(diff(xf))+abs(diff(yf)); I=find(d<=eps);
for i=length(I):-1:1, xf(I(i))=[]; yf(I(i))=[]; end

drawnow;

%% Compute lenght of line in cells
if nargin>3
    mode =2;
end
%% Right and left shore separately
% In this case, only matrix A is generated, holding the area of the stream
% as it exactly intersects with each cell
if nargin==3  && size(xyw,2)>2,  % Stream with separate left and right shore
    dxf=diff(xf);     dyf=diff(yf);     Lf =sqrt(dxf.^2+dyf.^2);
    E=[dxf./Lf dyf./Lf];  % direction vectors of all tracts
    e=[E(1,:);E(1:end-1,:)+E(2:end,:);E(end,:)];  % direction vectors at all corners, of length local w
    L=sqrt(e(:,1).^2+e(:,2).^2);
    e=[wf.*e(:,1)./L wf.*e(:,2)./L]/2;
     
    Rshore=[xf yf]+[+e(:,2), -e(:,1)]; % clockwise         = right shore points
    Lshore=[xf yf]+[-e(:,2), +e(:,1)]; % counterclockwise  = left  shore points
    % make right and leftshore one contour
    xf=[Rshore(:,1);Lshore(end:-1:1,1);Rshore(1,1)];
    yf=[Rshore(:,2);Lshore(end:-1:1,2);Rshore(1,2)];
    wf=zeros(size(xf));
    mode=3;
end
%% Intersect grid with polyline
xfm=([]); yfm=([]); wfm=([]);

% ==== intersect the interface with the grid
Np=length(xf);
if Np>1
    for ip=1:Np-1  % all vectors in interface points
        x1=xf(ip); x2=xf(ip+1);
        y1=yf(ip); y2=yf(ip+1);
        w1=wf(ip); w2=wf(ip+1);
        
        lam=unique([0; (x(:)-x1)./(x2-x1); (y-y1)./(y2-y1)]); I=find(lam<1 & lam>=0);
        xy=[x1+lam(I)*(x2-x1) y1+lam(I)*(y2-y1) w1+lam(I)*(w2-w1)];
        
        xfm=[xfm; xy(:,1)];
        yfm=[yfm; xy(:,2)];
        wfm=[wfm; xy(:,3)];
    end
    xfm=[xfm; xf(end)];
    yfm=[yfm; yf(end)];
    wfm=[wfm; wf(end)];
end
dxf=diff(xfm);
dyf=diff(yfm);

xm=0.5*(xfm(1:end-1)+xfm(2:end));
ym=0.5*(yfm(1:end-1)+yfm(2:end));
wm=0.5*(wfm(1:end-1)+wfm(2:end));

A=zeros(Ny,Nx);
for ip=1:length(xfm)-1
    ix=find(xm(ip)>x,1,'last');
    iy=find(ym(ip)<y,1,'last');
    
    switch mode
        case {1,3}
            if ~isempty(ix) && ix>=1 && ix<=Nx
                A(:,ix)=A(:,ix)+min(dy,max(0,y(1:end-1)-ym(ip))).*dxf(ip);
            end
        case 2
            if ~isempty(ix) && ~isempty(iy) && iy>=1 && iy<=Ny
                A(iy,ix)=A(iy,ix)+sqrt(dxf(ip).^2+dyf(ip).^2)*wm(ip);
            end
        otherwise
            error('Only case 1,2 or 3 allowed');
    end
end
A(abs(A)<1e-8)=0;


function selftest
x=0:25:1000;
y=(1000:-25:0)';

for i=1:length(y)
    plot(x([1,end]),y([i,i])); hold on;
end
for i=1:length(x)
    plot(x([i,i]),y([1,end]));
end
title('Selftest polyxgrid, mesh + lines');
drawnow

xyw=[
  118.6636  504.3860 30
  595.6221  214.9123 100
  570.2765  752.9240  30
  438.9401  656.4327 200
  289.1705  826.0234  20
];

[A,xf,yf]=polyxgrid(x,y,xyw);
plot(xf,yf,'-*r');
figure; spy(A); title('Selftest, shows cells that have been filled');

