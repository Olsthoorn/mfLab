% MFlow NWT - Newton a Newton Formulation for MODFLOW 2005
% Richard G. Niswonger, Sorab Panday, USGS 2011, Chapter 37 of section A
% groundwater Book 6: Modeling Techniques.
%
% Example Nr 2 from the manual
%
% TO 130425

clear variables; close all;

basename='MFlowNWT';

%% Example problem 2 peff=0.35;
% Simulation of a water-table mound, resulting from local recharge

%% variable/constants
kh     =    5;
kv     = 0.25;
Sy     = 0.20;
Ss     = 1e-5;
peff   = 0.35;
strthd =   25;
iPond  =    2;
d      = 5000/40;
Qpond  = 12500;
Apond  = 16*d^2;

t = [190, 708, 2630];

%% Grid
xGr=0:d:5000; 
yGr=0:d:5000;
zGr= 70:-5:0;

gr=gridObj(xGr,yGr,zGr);

%% Arrays
HK       = gr.const(kh);
VK       = gr.const(kv);
PEFF     = gr.const(peff);
SY       = gr.const(Sy);
SS       = gr.const(Ss);
STRTHD   = gr.const(strthd);
IBOUND   = gr.const(1); IBOUND(1:2,1:2,1)= iPond;
IBOUND(end,:,gr.zlay<strthd) = -1;
IBOUND(:,end,gr.zlay<strthd) = -1;

NPER = getPeriods(basename);
RECH = zeros(gr.Ny,gr.Nx); RECH(IBOUND(:,:,1)==iPond)=Qpond/Apond;

RECH = bsxfun(@times,ones(1,1,NPER),RECH);

save underneath peff
