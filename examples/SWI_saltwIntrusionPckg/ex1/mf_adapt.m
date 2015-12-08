%% Saltwater Intrusion Packge, example 1
%  Mark Bakker (2005)

%%
% A Rotating interface in a confined vertical cross section 40 m high and 250 m long
% and 2 m wide. There is a constant inflow of 2 m3/d from the left.
% The hydraulic parameteres are specified below.
% The interface is initially at 45 degrees angle from (x,z)=(80,0) to (x,z)=(120,-40).

%%
% There is one freshwater zone and one saltwater zone and thus one surface
% between them (NPNL=1). Flow is treated as stratified (ISTRAT=1).
% Since this is stratified flow with 2 zones, 2 density values must be specified
% for the SWI package in the MFLOW sheet (0 and 0.025).
% The source/sink terms are specified to be fresh water (zone 1); everywhere (ISOURCE=1);

basename='swiex1';

%% Specifics of this example
k      =    2;  % hydraulic conductivity
peff   =  0.2;  % effective porosity
strthd = 0.05;  % specified head in cell at right
Qo     =    2;  % source in cell 1

% for density sea SEAWAT worksheet
% for timestep see PER worksheet

%% Specify mesh
xGr=0:5:250;  % 50 cells 5 m wide
yGr=[0   2];  % widh is 1 m 
zGr=[0 -40];  % one layer

gr = gridObj(xGr,yGr,zGr); % request grid object

%% Model arrays

IBOUND  = gr.const(1);     IBOUND(:,end,1)=-1;  % fixed head at right hand side
HK      = gr.const(k);
VK      = gr.const(k);
STRTHD  = gr.const(strthd);  % start using fixed head
SSZ     = gr.const(peff);    % effective porosity called SSZ in SWI
ISOURCE = gr.const(2);       % sources are in zone 2 (saline, below the interface)

%% well
well = wellObj(basename,'wells',gr,HK,'PER'); % see sheet,' wells and Q1 in sheet PER for Q

%% ZETA planes
NSURF=1;  % We have only one surface here, so NSURF=1 in this case

x=[0 80 120 250];  % x points for interface
z=[0  0 -40 -40];  % z points for interface

ZETA = interp1(x,z,gr.xm);

