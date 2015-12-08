% MF_ADAPT -- Benchmark Problems and Application Example
% One-dimensional transport in a uniform flow field
% Zheng & Wang (1999) -- p130
%  TO 091106 091201


basename='OneD-Uniform';

%% Problem

% The problem considered here is transport due to uniform through a linear
% homogeneous model subject to dispersion, linear sorption and decay.
% Model has cells 10 m wide and is 101 cells long with fixed head at its
% right hand boundary and prescribed flow at the left end of the model.
% As the orginal problem only states the seepage velocity, prescribed flow
% at the left hand side of the model is the natural boundary.
% For this problem, conductivity is immaterial, any value is ok.
% We have the following cases
%
% aL= 0 m, R=0, lambda (decay) = 0
% aL=10 m, R=0, lambda (decay) = 0
% aL=10 m, R=5, lambda (decay) = 0
% aL=10 m, R=5, lambda (decay) = 0
%
% We simulate the four cases in a single model run by means of 4 layers
% that are separated by making sure that vertical dispersion is zero for
% each layer.
% Linear sorption is switched on with by setting ISOTHM=1 in worksheet
% MT3D. Decay is switched on by setting IREACT=1 in the same worksheet.
% The parameters to be set in the LAY worksheet are
%
%DSP    DSP    DSP    DSP      RCT   RCT    RCT    RCT    RCT
% AL    TRPT   TRPV   DMCOEF_1 RHOB  SP1_1  SP2_1  RC1_1  RC2_1
%  0     0.1      0      0        1      0      0      0      0
% 10     0.1      0      0        1      0      0      0      0
% 10     0.1      0      0        1      1      0      0      0
% 10     0.1      0      0        1      1      0  0.002      0
%
%
% The other parameters are defined below

%% Conductivities and effective porosity
kh= 10; kv=0;   % m/d
peff=0.25;      % []
vx=0.24;        % m/d,  % given seepage velocity = vx/peff=
q=vx*peff;      % Dary flux or specific discharge
C0=1;           % left boundary fixed concentration
% duration time - 2000 d, see worksheet PER

%% Mesh
NROW=1; NCOL=101; NLAY=4;

dx=10; dy=1; dz=1;

xGr=  dx*(0:NCOL);
yGr=  dy*(0:NROW);
zGr= -dz*(0:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% and to generate all other matrices
STCONC{1}=gr.const(0);

IBOUND=gr.const(1);  IBOUND(:,end,:)= -1;
ICBUND=gr.const(1);  ICBUND(:,end,:)= -1;
STRTHD=gr.const(0);

HK    =gr.const(kh);
VK    =gr.const(kv);
PEFF  =gr.const(peff);

%% East and west boundary cell index triples + period number (only 1 period)
iwest = 2;
ieast =-3;

IBOUND(:,  1,:)=iwest;
IBOUND(:,end,:)=ieast;

[WEL,PNTSRC]= bcnZone(basename,'WEL',IBOUND,[iwest,q],C0);

save underneath C0 peff