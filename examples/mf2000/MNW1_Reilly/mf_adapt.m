% MNW1 benchmark Halford and Hansen 2002 USGS open-file report 02-293
% The analyze part is still under construction and also the checking has
% yet to be completed TO 110822

%  TO 110808 110822

basename='MNW1_Reilly';

%% Mesh using table data

% Two aquifers separate by 50 ft confining bed.
% Top unconfined base 50ft, k=60 ft/d, Sy=0.05;;
% elevation bottom = 05ft. T2=15e3 ft2/d. lower quifer, S=0.001;
% 66 mi^2 21 rows and 14 columns of 2500 by 2500 ft cells.
% Specified head and drains in layer one and constant.

zGr=0:-5:-205;
yGr=sinespace(0,100,31,pi/5,pi/2); diff(yGr)
xGr=[sinespace(-250,-50,20,pi/2,pi/15) -47.5:2.5:47.5 sinespace(50,9750,200,pi/100,pi/2)]; diff(xGr)

gr = gridObj(xGr,yGr,zGr);

HK     = gr.const(250);   % ft/d full array neede as argument to mf_setmnwells
VK     = gr.const( 50);     % ft^2/d
STRTHD = gr.const(  0);     % ft
PEFF   = gr.const(0.2);
SY     = gr.const(0.05);       % -
SS     = gr.const(0.0001);     % -
%% Generate all other matrices using table data

IBOUND = gr.const(1); IBOUND(:,end,:)=-1;%  IBOUND(:,[1 end],:)=-1; IBOUND(1,:,:)=-1;

ICBUND =gr.const(01);
ICBUND(gr.ym<=50,gr.xm>=-25 & gr.xm<50,:)=1;

STCONC = gr.const(0);
STCONC(gr.ym<20,gr.xm>-20 & gr.xm<20,1)=100;  % contaminant

MNW = MNW1Obj(basename,'MNW',gr,HK,'Q_','C_');

save Underneath basename