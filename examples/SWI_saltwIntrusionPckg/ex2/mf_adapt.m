%% Saltwater Intrusion Packge, example 2

% The second examples is a variation of the first example, but includes three
% zones and no constant flow.
% The groundwater is divided into 3 zones: fresh, brackish, and salt water,
% with dimensionless desnities of 0, 0.0125 and 0.025, respectively.
% At first, flow is treated as stratified.
% Initially, at time t=0, both interfaces are straight and make a 45 degree
% angle with the horizontal and run from (x,z)=(150,0) to (x,z)=(190,-40),
% and from (x,z)=(110,0) to x,z)=(150,-40). The head is fixed to 0.05 at
% the top ceell 1. The brackish zonee will rotate to a horizontal position
% through time.

%%
% The Domain is discretized intno 60 cells of 5 m wide, and the timestep is
% 2 d. There are three zones and thus the initial elevations of two
% surfaces area specified.

%%
% The maximum slop of the toe and tip is chosen to be 0.4, and the
% deltazeta and deltazetamin parameters are computed according to equation 16.
%%
% All sources and sinks are fresh water, so ISOURCE==1.

%%
% For illustration puposes, the same porblem is solved with the variable
% density option. This requires two changes in the SWI input file:
% 1) the ISTRAT parameters must be set to 0 to indicated variable density flow
% (rather than stratified flow)
% 2) Second, the NU variable now needs 4 input values: the value of nu along
% the top of zone 1, the two surfaces and the for the bottom of zone 3, so that
% the fourth line in the swiex2.swi files becomes: 0 0 0.025 0.025.
% Notice that this gives an average value
% NU=0.0125 in the brackish zone.

basename='swiex2';

%% Specifics of this example
k      =    2;  % hydraulic conductivity
peff   =  0.2;  % effective porosity
strthd = 0.05;  % specified head in cell at right
Qo     =    2;  % source in cell 1

%% Specify mesh
xGr=0:5:300;
yGr=[0   2];
zGr=[0 -40];

gr = gridObj(xGr,yGr,zGr);

%% Matrices of computation grid
IBOUND  = gr.const(1); IBOUND(1,1,1)=-1;
HK      = gr.const(k);
VK      = gr.const(k);
STRTHD  = gr.const(strthd);
SSZ     = gr.const(peff);   % effective porosity called SSZ in SWI
ISOURCE = gr.const(1);

%% ZETA planes
surfaces = [interp1([0 150 190 300],[0  0  -40 -40],gr.xm); ...
            interp1([0 110 150 300],[0  0  -40 -40],gr.xm)];

ZETA = XS(surfaces);
