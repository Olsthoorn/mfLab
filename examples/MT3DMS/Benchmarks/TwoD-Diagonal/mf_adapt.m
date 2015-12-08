% MT3DMS Benchmark Problem see Zheng & Wang (1999), p139ff
% Two-dimensional Transport in a Diagonal Flow Field
%
%  TO 091114 091203 120416

basename = 'TwoD-Diagonal';
%GREP = 'STRESS';
RUNBACKGROUND = 1;

%% Model parameters specified by Zheng

[Pnms,Pvals] = getExcelData(basename,'model','Vertical');

NCOL = Pvals(strmatchi('NCOL',Pnms),1);
NROW = Pvals(strmatchi('NROW',Pnms),1);
NLAY = Pvals(strmatchi('NLAY',Pnms),1);
dx   = Pvals(strmatchi('DX'  ,Pnms),1);
dy   = Pvals(strmatchi('DY'  ,Pnms),1);
dz   = Pvals(strmatchi('DZ'  ,Pnms),1);
peff = Pvals(strmatchi('peff',Pnms),1);
k    = Pvals(strmatchi('k'   ,Pnms),1);
Q    = Pvals(strmatchi('Q'   ,Pnms),1);
cIn  = Pvals(strmatchi('cIn' ,Pnms),1);
cGrw = Pvals(strmatchi('cGrw',Pnms),1);
dhdx = Pvals(strmatchi('dhdx',Pnms),1);
dhdy = Pvals(strmatchi('dhdy',Pnms),1);

h0 = 0 ; % dummy head NW corner of model

%% Mesh

xGr= dx*(0:NCOL);
yGr= dy*(0:NROW);
zGr=-dz*(0:NLAY);

gr=gridObj(xGr,yGr,zGr);

%% arrays

STRTHD = h0+dhdx*(gr.XMlay-gr.xm(1))  + ...
           +dhdy*(gr.YMlay-gr.ym(1));   % gradient

IBOUND    = gr.const(1);
ICBUND    = gr.const(1);
STCONC{1} = gr.const(cGrw);  % Initial concentration of species 1
HK        = gr.const(k);
VK        = gr.const(k);  % dummy for one layer
PEFF      = gr.const(peff);

%% Boundary conditions
IBOUND(:,[1 end],:) = -1;
IBOUND([1 end],:,:) = -1;

iCCC = 4;   IBOUND(90,10,1)=iCCC;  % constand concentration cell

[~,PNTSRC] = bcnZone(basename,'CCC',IBOUND,{iCCC, Q},cIn);
