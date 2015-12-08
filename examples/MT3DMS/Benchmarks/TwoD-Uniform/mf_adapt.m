% MF_ADAPT --  mt3dms injeciton i n uniform flow field, example
% Zie Zheng & Wang (1999) p138ff
%
%  TO 091112

basename='TwoD-Uniform';

GREP = 'STRESS';

% Zheng (1999) benchmark example, page 138-139
% The problem deals with an injection well in a uniform flow field whose
% direction is at 45 degrees in NE direction. The analytical solution was
% given by Wilson and Miller (1978). The analytical solution is applicable
% under the condtion that the aquifer is of infinite exent and relatively
% thin in vertical direction, so that instantaneous vertical mixing can be
% assumed and the injection rate is insignificant with respect to the
% ambient uniform flow.

% The numerical model consists of 31 rows and 46 columns and 1 layer. The
% direction of the uniform ambient flow is along the rows of the model
% (x-direction). Teh model parameters are

%% Model parameters specified by Zheng

RUNBACKGROUND = 1;

%% Model parameters specified by Zheng

[Pnms,Pvals] = getExcelData(basename,'model','Vertical');

NCOL = Pvals(strmatchi('NCOL',Pnms,'exact'),1);
NROW = Pvals(strmatchi('NROW',Pnms,'exact'),1);
NLAY = Pvals(strmatchi('NLAY',Pnms,'exact'),1);
dx   = Pvals(strmatchi('DX'  ,Pnms,'exact'),1);
dy   = Pvals(strmatchi('DY'  ,Pnms,'exact'),1);
dz   = Pvals(strmatchi('DZ'  ,Pnms,'exact'),1);
peff = Pvals(strmatchi('peff',Pnms,'exact'),1);
vx   = Pvals(strmatchi('vx'  ,Pnms,'exact'),1);
vy   = Pvals(strmatchi('vx'  ,Pnms,'exact'),1);
k    = Pvals(strmatchi('k'   ,Pnms,'exact'),1);
Q    = Pvals(strmatchi('Q'   ,Pnms,'exact'),1);
cIn  = Pvals(strmatchi('cIn' ,Pnms,'exact'),1);
cGrw = Pvals(strmatchi('cGrw',Pnms,'exact'),1);

h0   =  10; % dummy head at west of model
qx   = vx*peff;
qy   = vy*peff;

%% Mesh definition
xGr= dx*(0:1:NCOL);
yGr= dy*(0:1:NROW);
zGr= dy*(0:1:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% Generate all other matrices
IBOUND    = gr.const(1);  IBOUND(:,end,:)= -1;  % fix right-hand head
ICBUND    = gr.const(1);   % no fixed concentration points needed

STRTHD = h0-qx*(gr.XMlay-gr.xm(1))/k;   % gradient

STCONC{1} = gr.const(cGrw);  % Initial concentration of species 1
STCONC{2} = gr.const(cGrw);  % Initial concentration of species 2

HK        = gr.const(k);
VK        = gr.const(k);  % dummy for one layer
PEFF      = gr.const(peff);

%% Input of boundaries through zones
[~,~,NPER] = getPeriods(basename);

iEast = 2; % east
iWest = 3; % west
iWell = 4; % well

IBOUND(:,1  ,:) = iWest;  % west fixed head
IBOUND(:,end,:) = iEast;  % east fixed head
IBOUND(16,16,1)=  iWell;  % well

zoneVals = { iWest, STRTHD(:,  1,:), STRTHD(:,  1,:);
             iEast, STRTHD(:,end,:), STRTHD(:,end,:)};
        
[WEL,PNTSRC              ] = bcnZone(basename,'WEL',IBOUND,{iWell, Q},cIn);
[CHD,PNTSRC(end+(1:NPER))] = bcnZone(basename,'CHD',IBOUND,zoneVals, cGrw);
