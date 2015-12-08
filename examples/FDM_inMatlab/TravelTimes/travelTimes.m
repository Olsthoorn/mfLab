% Set up a simple cross section model to compute the travel time
% TO 121115

%% Model values
L   = 350;   % widht of model
H   =  50;   % height of model
kh  = 10;    % horizontal k
por = 0.35;  % effective porosity
Q0  = -100;  % extraction [m3/d]
ddx = 2;     % cell size x-direction
ddy = 2;     % cell size z-direction
DZ  = 1;     % thickness of cross section, generally DZ=1 m
Eps = 0.001; % just a small number ot make very thin cells on all sides
             % to put the boundarys perfectly on the outside of the model
             % This is not necessary per se

%% Mesh
xGr = [-L/2-Eps -L/2:ddx:L/2  L/2+Eps];
yGr = [-H/2-Eps -H/2:ddy:H/2  H/2+Eps];

%% Cleand-up mesh
[xGr,yGr,xm,ym,dx,dy,Nx,Ny] = modelsize(xGr,yGr);

%% Model arrays
IBOUND = ones(Ny,Nx);       % nowadays like in MODFLOW (>0 computed <0 fixed, 0 inactive)
K   = kh  * ones(Ny,Nx);    % full K matrix
POR = por * ones(Ny,Nx);    % full POR matrix
FH  = 0   * ones(Ny,Nx);    % full starting or initial head matrix (used where fixed)
FQ  = 0   * ones(Ny,Nx);    % full fixed injection flows (0 if no injection or extraction)

IBOUND(:,[1 end])=-1;       % fix heads at left and right boundary of model 

% Set injection in center of model equal to Q0 (is negative, therefore extraction).
FQ  (ym>-ddy/2 & ym < ddy/2,xm >-ddx/2 & xm < ddx/2) = Q0;

%% Sun steady state 2D groundwater model (finite differences)
% Yields heads, net injection of all nodes, stream function at cell
% corners, right flow across cell faces and down flow acros cell faces
[Phi,Q,Psi,Qx,Qy] =fdm2(xGr,yGr,K,K,IBOUND,FH,FQ);

%% Show the heads and stream lines (stream lines are not flow paths, they may however conincide)
figure; hold on; xlabel('x [m]'); ylabel('z [m]'); set(gca,'xlim',xGr([2 end-1]),'ylim',yGr([end-1 2]));
title('Flow and travel time in fault (If Technology, greatings Theo Olsthoorn 15 Nov 2012')

prange = Q0*(-0.5:0.05:0.5);  % stream function values, I know the range: Q0

% Contour heads (just 25 lines, because I don't know the range beforehand
contour(xm,ym,Phi,25,'r');

% Plot the streamlines by contouring the stream function, using values in prange
contour(xGr(2:end-1),yGr,Psi,prange,'b');

% use same scales in z and x direction so that both are perpendicular
axis('equal');
axis('tight');

%% Travel times

% Starting points of flow paths
xSt = [ones(Ny-2,1) * xGr(2); ones(Ny-2,1) * xGr(end-1)];
ySt = [ym(2:end-1); ym(2:end-1)];

% times you want plotted
T = 10:10:200;

% Compute and plot path lines on existing figure. They will be marked using
% the marker ('.' in this case) and will be computed using times T.
% DZ is thickness of cross section, i.e. generally DZ=1.
% xGr,yGr is grid
% Q,Qx,Qy from fdm2 which was ran before.
[XP YP TP]=fdm2path(xGr,yGr,DZ,Q,Qx,Qy,por,T,'.',xSt,ySt);

% Each series of points of one flow path in XP,YP and YP is preceded and
% post-ceded with a NaN to easily find them
I=find(isnan(XP));          % so where are these NaNs ?
Ist  = I(1:end-1)+1;        % Then here start each path line
Iend = I(2:end  )-1;        % and here it ends
Nr   = (1:length(Ist))';    % generate path line numbers

% Plot total travel time
fprintf('Total travel times:\n');

% Header
fprintf('%-8s %-8s %-8s %-12s\n','Nr','xStart','zStart','time [d]');

% The total travel times
fprintf('%6d %8g %8g %12g\n',[Nr,XP(Ist),YP(Ist),TP(Iend)]');
