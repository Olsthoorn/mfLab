%% Tutorial for using the gridObj
% Teaches and demonstrates the use of the gridObj to generate and
% manipulate grid and 3D data arrays
%
% TO 151208

%% Objects in general and gridObj objects in _mfLab_
% Objects are "things" (instantiations) of a class defining its behavior and the
% data it can hold. Each Object may be regarded a container holding its class-specific
% data while all objects of the same class share common behavior. This behavior
% may be influenced by the actual data carried by a specific object.
% Grid objects belonging to the class gridObj are examples. Each grid object defines a finite difference grid.
% is an example of a class. Each instance may hold a different grid but all
% share behavior or functionality as defined by the methods (functions) of
% the class. For example a grid can provide the volume of each cell, area
% of cells, but can also plot itself.
%
% The gridObj is a class instantiations of which define a finite difference
% model grid or mesh. The gridObj to _mfLab_.
% _mfLab_ tries to make
% modelling as independent of a specific grid as feasible. Users should have
% to deal with grid specifics as little as possible.
% 
% To best practice with mfLab is to learn working with the grid object, observe
% its functionality and its use. Examples are given through this tutorial.
%
% Grid objects are always 3D as everything in _mfLab_.
% A finite difference grid is generated with the coordinates of the gridlines
% for the 3 directions specified: xGr, yGr, zGr. Extra parameters that may
% be used when generating a grid are LAYCBD and MINDZ.
% For now forget about LAYCBD and MINDZ.
% xGr, yGr and zGr define the coordinates of the grid lines, not the grid
% centers. The model is always aligned with its grid coordinates.

% Define the color for the background of the figures
grey = [0.8 0.8 0.8];

%% Lets get the size of our monitors
% Often we have more than one monitor, that we want to use to spread our
% figures on. This is useful when a large number of figures will be shown
% and the Matlab window is often covering the figure we already have
% plotted, while we wish to keep the visible to the extent possible.
% If we have more screens, the total combined screen area is called the
% canvas. It is continuous and the coordinates of the parts of it that are
% within each of our active displays are easily obtained with

scrSizes = get(0,'Monitorpositions');

% whicch for each screen gives the [xLL, yLL, width, height] values in pixels.
% With displays of different resolutions, some parts of the canvas cannot
% be seen. So one should stick within the limtis of the pixelcoordinates
% obtained for each display.

% lets get the default position of a generate figure
close all
figure; p = get(1,'position'); close(1);

%Now let's distribute NxM figures across one of the screens usig the
%default width and height 'wh' that Matlab just returned.

iDisp = 1; % use first display

N=5; % nr of figure columns
M=3; % nr of figure rows

H = round(scrSizes(iDisp,1) + scrSizes(iDisp,3) * (0:(N-1))/N);
V = round(scrSizes(iDisp,2) + scrSizes(iDisp,4) * (0:(M-1))/M);
pos = NaN(N*M,4);
k = 0;
for iv = M:-1:1
    for ih = 1:N
        k = k+1;
        pos(k,:) = [H(ih) V(iv) p(3:4)];
    end
end

% so we can put fiure k on position pos(k,:)
% figure('pos',pos(k,:));


%% Example of how to generate a gridObj
% Define gridlines in the three directions:
xGr = [-1000:100:2000 -250:25:250];
yGr = [-500:50:500 -200:20:200];
zGr = [ 0 -10 -20 -25 -40 -50 -80];

gr = gridObj(xGr,yGr,zGr); %#ok

% Notice that fully irregular grids may be obtained by specifying zGr as
% a full 3D array holding the tops and bottoms of all
% cells. These tops and bottoms can be different for each column,
% as is the case for a 3D Modflow model.

%% Example 3D grid of the Amsterdam Water Supply Groundwater model
% The 3D grid of the Amsterdam Water Supply Groundwater model is contained
% in the file data.mat which contains xGr,yGr and Z in model coorrdinates,
% i.e. local coordinates with the main y-axis parallel to the North-Sea
% coast.

%%
load data
gr = gridObj(xGr,yGr,Z);
%%
% Nothing seems to happen because the semi-colon suppresses the result
% being printed, but the grid object no exists.
% To show its contents one may type
display(gr);

%% 
% This shows the variable that the gridObj can provide a lot of information,
% all of which is based on the xGr, yGr and zGr (or Z array) used in the call.
%%
% Notice that in case a full grid Z is use in the call of the gridObj
% constructor, de vector zGr is actuall mean(mean(gr.Z),1,2), i.e a vector
% of size [1,1,Nz] with each value the layer-mean of the elevation of Z.
% So for plotting the actual grid, gr.Z must be used unless al layers are
% uniform. In that case gr.zGr and gr.Z are the same.

%%
% Show this grid as cross section along the x-axis (perpendicular to the
% North Sea coast along the row where yGr=0. This indeix of this row is
% selected using the function hit(yGr,0). See yGr for the available
% y-coordinates and try selecting other cross sections.
kFig = 1;
figure('pos',pos(kFig,:)); hold on;
title('x-Section through the grid of the Amsterdam Water Supply Model)'); xlabel('x'); ylabel('z');
gr.plotXSec(hit(yGr,0),1:gr.Nlay);
%%
% In the example above, the circumference of the cells were colored default
% using 'lines' 'on' as arguments. One may specify a different color such
% as 'llines' 'r' or 'lines' [1 0 1] as is standard in Matlab. One may also
% use different colors for the vertical lines and the horizontal lines by
% specifying 'hlines' and 'vlines' separately:
gr.plotXSec(1,1:gr.Nlay,'hlines','r','vlines','b');

%% The plotXSec has more options and can also be used in the y-direction
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
title('NS-Section through the grid of the Amsterdam Water Supply Model)'); xlabel('x'); ylabel('z');
gr.plotYSec(hit(xGr,0),1:gr.Nlay,'hlines',grey,'vlines',grey);

%% Instead of staircase lines around the cells as they are, continuous lines are possible
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
title('NS-Section through the grid of the Amsterdam Water Supply Model)'); xlabel('x'); ylabel('z');
gr.plotYSec(hit(xGr,0),1:gr.Nlay,'lines','k','smooth');

%%
% gr is an instance of the class gridObj. It hold concrete data of this
% grid, xGr, yGr and zGr, but shares behavior with any other instance of
% this class.
% This behavior is implied by the method (functions) defined for this class.
% At the bottom in blue you find the word "methods". Clicking on it yields a list
% of the methods defined for the class gridObj.
% To see the definition of any of the methods typ the grid class name
% followed by a dot and the method:
eval('help gridObj.plotgrid');

%%
% To inspect the code in detail, pull it into the editor:
edit gridObj.plotgrid;

%%
% Plotgrid is an example of a method of the class gridObj: It knows how to plot a grid
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
xlabel('x [m]'); ylabel('y [m]'); title('Grid of AWS groundwater model');
gr.plotGrid(gca,'iy',1);

    %%
% Grid can provide a lot of information, which is computed upon request.
% This may be simple numbers like:

%%
gr.Nlay % Number of layers in the grid
%%
gr.Nx   % Number of columns (x-direction)
%%
gr.Ny   % Number of rows (y-direction)
%%
gr.size % size of the grid (Ny,Nx,Nlay) in natural order of Matlab
%%
gr.xm   % cell center coordinates
%%
gr.dx   % cell widths

%%
% Capital letters are used for arrays as opposed to vectors. Compare
gr.XM;   % cell centers of all cells of the 3D grid
gr.Xm;   % cell centers of all cells in a layer (2D instead of 3D).
gr.AREA; % the area of all cells in a horizontal plane
%% Surface area of entire grid
gr.area
%% Volume of subgrid cells
% To obtain the volume of the individual cells in part of the grid specify the
% indices of the desired part. You may use matlab function between to
% easily extract a part by specifying coordinates instead of indices
gr.Vlay(5:10,between(gr.xm,1000,2000),2:5)

%% Grid line coordinates are sorted and made unique
% While the grid object is generated, the grid lines are sorted and made
% unique by eliminating repeated values. The coordinates in the z-direction
% will be sorted from high to low. The first layer will always correspond
% to the top of the model. One must ensure that the data for the
% y coordinates are sorted high to low as well. This ensures that the layout
% of printed layer of grid data has the same orientation as in a plan view of
% reality. However, the line with the lowest coordinate then corresponds
% with y(end) and not with y(1), unless we deal with cross sections
% consisting of a grid with a single row.
%
% There are some functions useful for specifying grid coordianates such as
% Matlab's linspace and logspace and _mfLab_'s sinespace. Sinespace divides a
% distance between two points, say a and b into n-1 cells, where the
% width of these cells varies according to the sine function. The initial
% and final width of the cells can be adjusted by specifying the start and
% finish sine by means of the start and finish angle angle1 and angle2 in
% radians.
%%
eval('help sinespace');

%%
% To divide the space between 0 and 100 into 50 cells where the width varies
% according to the sine function between pi/10 and pi/2:
xGr2 = sinespace(0, 100, 51, pi/10, pi/2); %#ok
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
title('Demonstrating the use of sinespace');
xlabel('x');
dy=0.25;

iDisp=0; x0 = 50;

propvals={'HorizontalAlignment','center'};
% Use sinespace for instance to generate a grid that starts and ends with
% small cell widths and has broad distances between points in the center.
% Or use it to start with almost constant cell width (first angle pi/2) and
% end fine (second angle pi).
xGr2 = sinespace(0, 100, 51, pi/2, pi);
iDisp=iDisp+1; plot(xGr2,iDisp*dy* ones(size(xGr2)),[mf_color(iDisp),'.']);
text(x0,(iDisp+0.25)*dy,'sinespace(0, 100, 51, pi/2, pi)',propvals{:});

xGr2 = sinespace(0, 100, 51, 0, pi);
iDisp=iDisp+1; plot(xGr2,iDisp*dy* ones(size(xGr2)),[mf_color(iDisp),'.']);
text(x0,(iDisp+0.25)*dy,'sinespace(0, 100, 51, 0, pi)',propvals{:});

xGr2 = sinespace(0, 100, 51, -pi/2, pi/2);
iDisp=iDisp+1; plot(xGr2,iDisp*dy* ones(size(xGr2)),[mf_color(iDisp),'.']);
text(x0,(iDisp+0.25)*dy,'sinespace(0, 100, 51, -pi/2, pi)',propvals{:});

xGr2 = sinespace(0, 100, 51, -pi, pi);
iDisp=iDisp+1;plot(xGr2,iDisp*dy* ones(size(xGr2)),[mf_color(iDisp),'.']);
text(x0,(iDisp+0.25)*dy,'sinespace(0, 100, 51, -pi, pi)',propvals{:});

iDisp=iDisp+1;
set(gca,'ylim',[0 iDisp]*dy);

%% Further grid-generating arguments LAYCBD, MINDZ and AXIAL.
% LAYCBD is a vector specifying whether a model layer (LAY) has a confining
% bed (CBD) below it. Only model layers have cells and 
% heads and flow terms computed in and for them. A confining bed is a
% "virtual" hydraulic resistant layer for which only the leakance or
% the vertical conductivity is specified, but without heads.
% CBD layers have no cells and do not count as model layers;
% they determine the hydraulic leakance (1/resistance) between tow model
% overlying modellayers.
% Yet, CBD layers do have a thickness and, therefore, their top and bottom
% elevation has to be specified.
% _mfLab_ requiers the elevation of all tops and bottoms
% to be specified in a single array zGr as input argument to the grid object.
% Confining beds represented the only method originally available
% in MODFLOW before MODFLOW 2000 (BCF-package). The leakance of the CBD layers
% was called VCONT [1/time]. VCONT is the reciprocal of the resistance c
% [time] of aquitards  that is often used outside the USA.
% The VCONT concept was used in quasi-3D modeling, where generally,each
% model layer represented an aquifer and each confining bed the aquitard
% between two overlying aquifers.
% This option is still present in the modern LPF package available since MODFLOW 2000.
% In the LPF package VKCB (vertical k of confining bed) are specified
% instead of the previously mentioned VCONT.
% In LPF package, confining beds may also be simulated without making use
% of the confining-bed concept, simply by using a model layer
% with low vertical conductivity. This method may be more intuitive and
% is necessary if information within confining beds is required. On the other
% hand, it adds extra layers with cells to be computed and thus requires
% more computation time. For transport modelling, like with MT3DMS or
% SEAWAT, it is advised to always represent confining beds by one or more model
% layers with low vertical conductivity, as the transport simulator cannot
% store and account correctly for the species in a confining bed.
%%
% The LAYCBD vector has one value for each model layer. If the layer in
% question has a confining bed attached to its bottom, the value is 1, and
% zero if not. This implies that the total number of layers and
% confining beds equals:
%%
display(length(gr.LAYCBD)+sum(gr.LAYCBD==1));
%%
% The first term equals the number of model layers and the second term
% the number of confining beds.
% We did not specify LAYCBD when generating the grid for the AWS
% groundwater model. Therefore it has all zeros
display(gr.LAYCBD);
%%
% This implies that all layers are model layers with computable cells.
% was not specified in the call of the gridObj constructor, it is
% considered to be all zeros, i.e. that there are no confining beds. To obtain the
% current LAYCBD vector, type
%%
% The actual LAYCBD used in the AWS model is
display(LAYCBD);
%%
% Each value represents a model (cell) layer, of which we have 10. Each
% value of 1 represents a cell layer + a CFB attached to its bottom. So we
% have 7 of those. The last value of LAYCBD is always zero.
% The Z array represents the tops of the model followd by the bottom of all
% layers including those of the confinging beds. So every zero in LAYCBD
% adds one laye in the Z-array and every 1 adds two layers in the Z-array.
display(length(LAYCBD)+sum(LAYCBD==1));
%%
% which adds up to 17
%%
% Without confining beds, Nlay == Nz, where Nlay is the number of layers
% (with cells) and Nz the number of layers + confining beds.
% Using the original grid again:
%%
display(gr.Nlay);
%%
display(gr.Nz);
%%
gr.Ncbd

%%
% Let us regenerate the original grid, but with the actual LAYCBD
gr = gridObj(xGr,yGr,Z,LAYCBD);
%%
% Notice that the number of layers is now smaller and that the number of
% confining beds is no long zero
%%
gr.Nlay
%%
gr.Nz
%%
gr.Ncbd

%%
% show a cross section:
y_section = 0;
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
set(gca,'xlim',gr.xc([1 end]),'ylim',[-250 50]);
gr.plotXSec(hit(gr.yGr,y_section));

%%
% It is also possible to show the grid in a transparent wire frame
% You can rotate the figure using the rotation tool in the menu bar of the
% matlab figure
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;

set(gca,'cameraposition', [-25000 -25000 300]);
set(gca,'cameratarget',   [     0      0 -150]);

gr.plotMesh; % wireframe mesh

%%
% Or we can plot the layers in combination with the mesh
kFig=kFig+1;
figure('pos',pos(kFig,:)); hold on;
xlabel('x'); ylabel('y'); zlabel('z'); title('layer and mesh ASW model, colors are layer elevations');

set(gca,'cameraposition', [-25000 -25000 300]);
set(gca,'cameratarget',   [     0      0 -150]);

%gr.plotLayers(1:gr.Nz,gr.ZM,'edgecolor','none');
gr.plotMesh;
colorbar
% Note that you can rotate in 3D using the tool in the menu of the
% resulting Matlab figure.

%% If you have wells of class wellObj you can plot them

%well = well.plotXY(gca);

%%
% Without having specified anything but the row along which the cross
% section is to be plotted, only the confining beds are shown; no other
% lines are drawn.
% To plot the same cross section with the horizontal
% interface between all layers use:
kFig=kFig+1;
figure('pos',pos(kFig,:)); hold on;
set(gca,'xlim',gr.xc([1 end]),'ylim',[-250 50]);
xlabel('x [m]'); ylabel('z [m MSL]');
title(sprintf('XS AWD model perpendicular to the coast at y=%.0f',y_section));
grey=get(gcf,'color');  % get background grey color from figure

gr.plotXSec(hit(gr.yGr,0),(1:gr.Nlay),'hlines',grey);
%%
% In the west are the contours visible of the young dunes (white) and in
% the east those of the old dunes (blue) which are represented by a
% different model layer.

%%
% To extract the full-size LAYCBD vector:
gr.LAYCBD
%%
% To replace LAYCBD with a new vector to reconfigure the choice of model layers
% versus confining beds, just regnerate the grid reusing the original
% grid input but now with the new LAYCBD. The LAYCBD vector needs not to be
% complete, it suffices to specify the vector down to the lowest layer with
% a confining bed. mfLab will complete it and ensures the last layer gets
% value zero.
%%
% To ensure that the internal proceures do not set LAYCBD directly like this:
gr.LAYCBD = LAYCBD; % strongly discouraged!
%%
% But regenerate the gridObj in its entirety
gr2 = gridObj(gr.xGr(50:100),gr.yGr(20:30),gr.Z(20:30-1,50:100-1,:),LAYCBD);

%% Generating grid arrays using method gridObj.const
% The grid object is very useful (and elegant) in generating 3D grid arrays.
% These arrays can be generated using the gridObj method const.
% For example, generating a 3D array IBOUND of the size of the grid (cells)
% with all values set to 1, do this:
IBOUND = gr2.const(1);
%%
% However if the argument is a vector, the 3D array will have a many layers
% as the length of the vector each layer having the corresponding value of
% the vector. So if our 4 layer and 2 cbd model has the following values
% per layer and confining bed:
%%
kh   = [10 5 25 20 30 25 20 10 25 20];         % 10 layers
vkcb = [0.01 0.02 0.01 0.025 0.02 0.03 0.015]; % 7 CFB

% We can generate full arrays of this layer horizontal conductivity and the
% confining bed vertical conductivity using
%%
HK   = gr2.const(kh);     % remove semi colon to show
XS(HK(1,:,:))            % show in cross section
%%
VKCB = gr2.const(vkcb);   % remove semi colon to show
XS(VKCB(1,:,:))

%% other useful parameters are tops and bottoms of layers and confining beds
gr2.ZTlay;  % remove semi colon to show
gr2.ZBlay;  % remove semi colon to show
gr2.ZTcbd;  % remove semi colon to show
gr2.ZBcbd;  % remove semi colon to show


%% Advanced use of plotXSec
% A more advance use of gr.plotXSec is to not only plot/color the confining
% beds % (if they exist in the model) but also layers that have a
% conductivity less then some value k_low say
k_low = 0.5;

kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on
title('Cross section through grid'); xlabel('x'); ylabel('z');
gr2.plotXSec(1,mean(mean(HK,1),2)>k_low,'lines',grey,'linewidth',2);

%%
% The second argument is now a logical vector in z-dimension that is true
% for leays fulfilling the condition. Note that mean(...,1) is taking the
% mean along the y-direction and mean(...,2) along the x-direction. So each
% value of this vector is the mean of a layer. See help mean.

% Many more options can be provided using propertyName propertyValue pairs
% as further arguments.

%% Switching between 2D and 3D
% The gridObj dealw with 3D only. There is no 2D version. Everything
% in _mfLab_ is 3D as is in fact the case with the underlying groundwater
% codes. To facilate switching between 2D and 3D that is, by interchanging
% the first (y) and third (z) dimension using function XS(...)
%%
XS(HK(1,:,:))
%%
XS(VKCB(1,:,:))
%%
% XS also works the other way around, i.e. switching from a 2D vertical
% plane(Nz,Nx) to a 3D plane(1,Nx,Nz). Its implementation is trivially simple:
% function A = XS(A), A=permute(A,[3,2,1]); end
% To switch between vertical cross sections along the y-axis and the 3D
% situation, use YS
%%
YS(HK(:,1,:))
%%
YS(VKCB(:,1,:))

%% MINDZ
% The grid constructor can be called with a 5th argument MINDZ, which is
% the minimum layer thickness in the model. This will then be guaranteed.
% To reques MINDZ
%%
gr2.MINDZ

%%
MINDZ=0.01;

% To set a the grid with a new MINDZ use
gr3 = gridObj(gr.xGr,gr.yGr,gr.zGr,'MINDZ',MINDZ);

% This call forces the gridObj to carry out any internal procedures
% necessary to guerantee the integrety of the resulting grid. This is
% opposed to
gr3.MINDZ= MINDZ;

%% AXIAL
% The last, 6th argument of the gridObj constructor is AXIAL. If axial is
% non-zero, the model will run in axial-symmetric mode without any further
% action being required by the user. _mfLab_ knows which arrays it has to multiply
% by 2piR to achieve this. Hence, a single number AXIAL=0 versus AXIAL=1 as
% last gridObj argument, makes for the difference between a flat
% cross section and an aial-symmetric one. This is probably a unique
% feature.
%
% To reques AXIAL
gr3.AXIAL

%%
% To set axial
gr3.AXIAL=1;
% You may further explore the grid object yourself.
% You may even add your own functionality to it.

%% Generating starting points for particle tracking with MODPATH
% The grid in combination with a zoneArray can provide starting point
% locations for MODPATH. The call is
% STARTLOC = gr.startLoc(zoneArray,zoneVals)
% where zoneVals is a cell array having on each line
% zonVals{iz,:) = {zoneArray, izone, n, iFac, t_released}
%  izone is zone number, which must exist in zoneArray
%  n is [Nx Ny Nz] or just n of Nx=Ny=Nz the distribution of starting
%     points along iface or over the cell in 3D
% iFace is [iFace1 iFace2 .... ] to simultaneously generate starting points
%     for several iFaces (see MODPATH manual). Use zero to distributed in
%     3D within the cells. Use [1 2] to distribute Left-Right [3 4] to
%     distributed [Front-Back] [5 6] for [Top-Bottom] or any combination
%     which could be just one iFface or all 6 [1 2 3 4 5 6].
%     Notice that if iFace is used, the points will be released at the
%     particular cell face, if iFace=0 the points will never be at a cell
%     face, they will always be released inside the cell.
%     so with t_releassed = 365, n=5 and iFace=2 in the cells beloning to the
%     right boundary of the model 
%%
iCircle1 = 66;  % zone 1 index
iCircle2 = 77;  % zone 2 index
%%
% Distance from center of model to all cells of model
Dist = gr2.dist(mean(gr2.xGr),mean(gr2.yGr));
%%
% put zone number iCircle1 in top layer of IBOUND
IB1 = IBOUND(:,:,1);                  % separate top layer of IBOUND for a moment
IB1(Dist(:,:,1)<500)=iCircle1;     % insert zone iCircle1 in this layer only
IBOUND(:,:,1)=IB1;                    % replace top layer of IBOUND with IB1

%%
% put zone number iCircle2 in bottom layer of IBOUND
IB2 = IBOUND(:,:,end);
IB2(Dist(:,:,end)<350)=iCircle2;
IBOUND(:,:,end)=IB2;

% Specify values for two zones iCircle1 and iCircle2, for example
n1     = 5;      % 5 points per iface
n2     = [2 3 5]; % points per face along x,y and z axes
iFace1 = 2;      %   East
iFace2 = [3 5];   % front + back
t_release1 = 365.25; % 1 year
t_release2 = 730.50; % 2 years

zoneVals = {iCircle1 n1 iFace1 t_release1;...
            iCircle2 n2 iFace2 t_release2};
        
LOCATIONS = gr2.startLoc(IBOUND,zoneVals);

modelLoc = gr2.relloc2model(LOCATIONS);
kFig = kFig+1;
figure('pos',pos(kFig,:)); hold on;
view(3); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
for izone=1:length(modelLoc)
    plot3(modelLoc{izone}(:,1),modelLoc{izone}(:,2),modelLoc{izone}(:,3),[mf_color(izone) 'o']);
end




%% Check grid generation using a simple grid  of 3 by 3 by 3

% Generate a small grid for the sake of clarity
gr4 = gridObj(0:4,0:4,0:4); % 4x4x4 cells

IBOUND = gr4.const(1);

% Choose a zone number
izone1=2;

% Put zone into IBOUND
IBOUND(2,2,2)=izone1;

% Specify values for start points using n and  iFace
% n determins the distribution of points across iFace
% for iFace see MODPATH manual, allow more than one iFace specified
% simultaneously
n1     = 3;       iFace1 = 1;    % means [Ny=Ny=Nz=1]     iFace1=West;
n2     = 2;       iFace2 = 2;    % means [Nx=Ny=Nz=2]     iFace2=East;
n3     = [3 0 4]; iFace3 = [3 4];% means [Nx=2 Ny=2 Nz=4] iFace3=[Front Back]
n4     = [2 3 3]; iFace4 = [5 6];% means [Nx=5 Ny=6 Nz=7] iFace4=[Bottom Top]

%% iFace = 0 implies distribution of points in the interior of the cell
n5     =  3;      iFace5 =  0;   % means [Nx=3 Ny=3 Nz=3] iFace=0 (distributed throughout cell)

t_release = 365.25; % 1 year

%% Put all these starting points in the same zone
zoneVals = {izone1 n1 iFace1 t_release;...
            izone1 n2 iFace2 t_release;...
            izone1 n3 iFace3 t_release;...
            izone1 n4 iFace4 t_release;...
            izone1 n5 iFace5 t_release};

%% Get their relative locations in the from as required by MODPATH        
LOCATIONS = gr4.startLoc(IBOUND,zoneVals);
display(LOCATIONS{end}); % show the last package of these locations

%% Plot the starting points in the grid
kFig=kFig+1;
figure('pos',pos(kFig,:)); hold on;
view(3); xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Small grid with starting points for MODPATH generated by different methods');
%%
% Plot a mesh of the grid
gr4.plotMesh;

%%
% compute model coordinates of the locations
modelLoc = gr4.relloc2model(LOCATIONS);

%%
% the starting points had been packaged in cells with the points of each
% call in its own cell array. For all these bundles, plot the points in
% them with a distinctive color
for izone=1:length(modelLoc)
    plot3(modelLoc{izone}(:,1),modelLoc{izone}(:,2),modelLoc{izone}(:,3),[mf_color(izone) 'o']);
end



%% Brendan O'Boyle
kFig = kFig+1;
figure('name','Brendan O''Boyle''s figure','pos',pos(kFig,:)); hold on;
title('x-Section through the grid of the Amsterdam Water Supply Model)'); xlabel('x'); ylabel('z');
gr.plotXSec(hit(yGr,0),1:gr.Nlay,'lines','on');

