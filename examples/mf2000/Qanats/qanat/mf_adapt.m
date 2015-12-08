% Qanat, water or drainage tunnels
%
% This example concerns a qanat. It is insprired by Spaks and Smith (1986)
% and Lightfoot (1996) discussion about the Khettara in the Tafilalt Region
% in Marocco. In the Tafilalt valley many of the old Khettara are easily
% recognized on Google earth. The Tafilalt Valley can be readily foudn in
% google earth be looking for Tafilalt, Maroc. Doing so, one recognizes the
% map drawn by Dale Lightfoot in 1996..
%   --- climate
%   --- leakage of a qanat tunnel
%   --- yearly variation and drought periods response
%   --- effect of growth of nearby tube well numbers
%   --- effects of maintenance
%   --- effects of neighboring qanats
% We may also use the model to investiate whch rules should be installed
% to guarantee sustainability.
% The estimated flow can be confronted with demand and the size of the
% population that can be sustained.
%
% Our model is inspired on Parks & Smith (1983) Factors affecting the flow
% of aflaj in Oman: A modeling study. Journal of Hydrology, Vol 65, No. 4.
%
% Their modeling can now readily be carried out with MODFLOW. Yet we have
% to be very suspicious with respect to how MODFLOW handles dry cells. The
% MODFLOW way is a disaster. Dry cells are taken out of the model and also
% no discharge from a dry-cell area is computed, meaning that areas where
% dry cells occur in MODFLOW get no recharge!! This is definitely wrong.
% Therefore, we prefer using John Doherty's MF2KASP version of modflow,
% which prevents cells from running dry altogether, and, therefore, no
% water is actually lost.

clear variables; close all
basename='qanat';

%% Model grid
% The grid should be rpresentative of a mountanous foreland, which collects
% recharge and is suitable for constructing the qanat. Hence, we consider a
% water-table aquifer bounded at a water divide (closed) boundary on all
% sides, with a drainage area near the discharge side to incorporagte
% groundwater outflow. This boundary may be ground evaporation or marshland
% or a stream.

%% Basic parameters of the model and its grid

% The qanat will be oriented along the x-axis. The yAxis is considered
% parallel to the mountain range. The outlet of the qanat is located 
% at [x,y,z] [0,0,0]; The mountain range is at x=LX=8000; North and South
% boundaries are at LY and -LY respectively.
% The aqufier is considered of zero thikness at x=LX, increasing in
% thickness to 100 m at x=0. This reflext a slowly declining bedrock.
% Ground surface at x=LX is +40. The inclination of the 6 km long qanat is
% 1:1000, hence its toe coordinates are x=6000, z=6 m.

%% Improvements
%
% Downstream boundary should not be fixed
% Needs transient flow to deal with dry rivers in summer
% Wetting capabilities/dry cell with MODFLOW is kind of hopeless
% Perhaps the new NWT package is better.
% Dohrties ASP package can be used under windows. that works a million
% times better dan USGS's MODFLOW for like situations.
% TO 130323
% 
Nz = 5;

LX  =8000;  % [m] lengh of aquifer upstream of Qanat outlet
LY  =2000;  % [m] width of aquifer perpendicular of Qanat (wadi width)

zT  =   0;   % [m] elevation of ground surface at outlet of qanat
zB  =-100;   % [m] elevation of bedrock at outlet of qanat
zEnd=  40;   % [m] elevation of bedrock outcrop

dx  = 200;   % [m] cell size in x-direction
dy  =  50;   % [m] cell size in y-direction

k    =10;    % [m/d] conductivity
peff=0.35;   % [-]   effective porosity
cDrn = 1000;  % [m2/d] drain conductance
sy  =0.25;   % [-]   specific yield
ss  =1e-5;   % [1/m] specific elastic storage

%% The model grid
xGr   = 0:dx:LX;      xm = 0.5*(xGr(1:end-1)+xGr(2:end));
yGr   = LY:-dy:-LY;   ym = 0.5*(yGr(1:end-1)+yGr(2:end));

zT    = ones(size(ym(:)))*interp1(xGr([1 end]),[0  zEnd],xm);
zB    = ones(size(ym(:)))*interp1(xGr([1 end]),[zB zEnd],xm);
DZ    = (zT-zB)/Nz;

Z = NaN(numel(ym),numel(xm),Nz);
Z(:,:,1) = zT;
for iz=1:Nz, Z(:,:,iz+1) = Z(:,:,iz)-DZ; end

gr = gridObj(xGr,yGr,Z);

%% IBOUND array all ones as there are no fixed head boundaries
IBOUND =  gr.const(1);  % all ones, i.e. ordinary cells
IBOUND(:,1,:)=-1;

STRTHD =  gr.const(50);  % sufficiently hight to get a good start
STRTHD(:,1,:) = -5;
%% HK and VK
HK = gr.const(k);
VK = gr.const(k/5);

%qanat = drnObj(basneame,'drains',gr.HK);
zQanat = interp1([0 6000],[0 6],gr.xm);
xQanat = gr.xm(~isnan(zQanat));
zQanat = zQanat(~isnan(zQanat));
yQanat = zeros(size(xQanat));
Iper   = ones(size(zQanat(:)));
cQanat = cDrn * ones(size(zQanat));
[Ix,Iy,Iz] = xyzindex(xQanat(:),yQanat(:),zQanat(:),gr);

DRN = [Iper,Iz,Iy,Ix,zQanat(:),cQanat(:)];

DRNidx = cellIndex(DRN(:,[4,3,2]),gr);

save underneath DRNidx xQanat yQanat zQanat  % needed in mf_analyze to extract observations
