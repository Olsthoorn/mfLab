% Cross section showing functioning Large-Scale storage of desalinated
% water in Abu Dhabi. Photos are from field excursion in October 2010
% with the ISMAR conference.
%
% TO 071024 071127 110316 100327 110413
%
clear variables;
close all;

AXIAL=0;

basename='Oasen';

peff=0.35;
ss  =1e-5;
sy = 0.1;
cfresh=50/20000;  % relative salinity concentration
RivBotElev=-5; % [m NAP] elevation of river bottom
cRiver  = 30;  % [d] resistance of river bottom
cPolder = 50;  % [d] resistance of cover layer (eeklaag polder)
hPolder = -1.5; % maintained water level in the polder

%% Get configuration from excel workbook

xGr=-1500:10:1500;           % x grid lines [m]
yGr=[-0.5 0.5];              % y grid lines [m] make it a cross section of 1 m width;
zGr=[-1.5 -13 -23 -50   0:-5:-50 -50:-10:-100];   % add sufficient lines to make a sufficiently fine grid

[xGr,yGr,zGr,xm,ym,zm,Dx,Dy,Dz,Nx,Ny,Nz]=modelsize3(xGr,yGr,zGr); % coordinates housekeeping

% generate full 3D matrices of cell sizes and cell center coordinates for layer use
[DX,DY,DZ]=meshgrid(Dx,Dy,Dz);
[XM,YM,ZM]=meshgrid(xm,ym,zm);

Z=repmat(zGr,[Ny,Nx,1]); % full Z array

 HK       =mf_zone(basename,xGr,yGr,zGr,'Config','Material','kh'); % horizontal conductivities
 [VK,Conf]=mf_zone(basename,xGr,yGr,zGr,'Config','Material','kv'); % vertical   conductivities

%% fixed head boundaries. Are interpreted as local point heads in fdm2dens
IBOUND=ones(Ny,Nx,Nz);

IBOUND(:,[1 end],:)=-1;  % left and right hand side are fixed head boundaries


%% We will use IBOUND to mark the cells that will be GHB and RIV
%  GHB for polder cells, representing cover layer resistance
%  RIV for river  cells

IBOUND(:,:,1)=3;  % use 3 for GHB

%% We will use RIV to fix river stage (river levels)

xRivL=Conf.xL(strmatchi('River',Conf.names)); % x of left  side of river
xRivR=Conf.xR(strmatchi('River',Conf.names)); % x of right side of river

% see that we overwrite IBOUND==3 where the river is
IBOUND(:,xm>xRivL & xm<xRivR,1)=2;   % use 2 to indicate river cells in first layer

%find(IBOUND==2)  % check that IBOUND==2 exists (== riv cells)

PEFF = ones(Ny,Nx,Nz)*peff;
SS  =  ones(Ny,Nx,Nz)*ss;
SY  =  ones(Ny,Nx,Nz)*sy;

%% Start heads and concentrations

STCONC=ones(Ny,Nx,Nz)*cfresh;
STRTHD=zeros(Ny,Nx,Nz);

%% Extractions. 

[wel,WEL]=mf_setwells(basename,xGr,yGr,Z,HK,'wells');



%% boundary conditions

Ighb=find(IBOUND==3);  % global indices of river cells
Iriv=find(IBOUND==2);  % global indices of river cells

Ighb=Ighb(:);
Iriv=Iriv(:);

LRCghb=cellIndices(Ighb,size(IBOUND),'LRC'); % Layer Row Col of GHB points
LRCriv=cellIndices(Iriv,size(IBOUND),'LRC'); % Layer Row Col of RIV points

ughb=ones(size(Ighb)); % column vector of ones; length=length GHB points
uriv=ones(size(Iriv)); % column vector of ones; length=length RIV points

ghbCond = DX(Ighb).*DY(Ighb)/cPolder;
rivCond = DX(Iriv).*DY(Iriv)/cRiver;  % compute river bottom conductance

[pernams,pervals]=getPeriods(basename); % read stress period information (need rivers stage)

RivStage=pervals(:,strmatchi('RivStage',pernams));
NPER=size(pervals,1);

GHB=[];
RIV=[];
for iPer=1:NPER
    GHB=[GHB; ...
        ughb*iPer LRCghb ughb*hPolder        ghbCond];
    RIV=[RIV; ...
        uriv*iPer LRCriv uriv*RivStage(iPer) rivCond uriv*RivBotElev];
end

%% River package

% zie boven. De rivier is nu ingebouwd via IBOUND==2 als rivier cellen


% Nriv=20;
% 
% Nrowriv=ones(Nriv,1);
% for i=2:Ny;
% Nrowrivnew=ones(Nriv,1)*i;
% Nrowriv=[Nrowriv; Nrowrivnew];
% end
% 
% Ncolriv=1:Nriv;
% Ncolriv=repmat(Ncolriv',length(Nrowriv)/Nriv,1);
% 
% stage=ones(length(Nrowriv),1,1)*0;
% cond=ones(length(Nrowriv),1,1)*(25/0.5)*Dx*Dy;
% Rbot=ones(length(Nrowriv),1,1)*-12;
% 
% RIV=[ones(Ny*Nriv,2) Nrowriv Ncolriv stage cond Rbot];
% 
% for iPer =2:NPER
%     Rivnew=[ones(length(Nrowriv),1)*iPer ones(length(Nrowriv),1) Nrowriv Ncolriv stage cond Rbot];
%     RIV=[RIV; Rivnew];
% end

%% Point source in river
% We have to define point sources, that is prescribe the concentration of water entering the model
% through all boundary points. That is if this concentration is nonzero. If it is zero, then we may
% skip this.

% all concentrations are taken relative to the csalt, this is why
% the concentrations are multiplied by csalt

ItypeWel=2; % see MT3DMS manual (p122)
ItypeGhb=5; % see MT3DMS manual (p122)
ItypeRiv=4; % see MT3DMS manual (p122)
concRiv = pervals(:,strmatchi('cRiv',pernams)); % conc river water
concPol = 0;  % concentration polder water
PNTSRC=[];
for iPer=1:NPER
    for iw=1:length(wel)
        uw=ones(size(wel(iw).idx));
        PNTSRC=[PNTSRC; ...
            uw(:)*iPer, wel(iw).LRC ...
            uw(:)*pervals(iPer,strmatchi(['C1_' sprintf('%d',iw)],pernams)), ...
            uw(:)*ItypeWel];
    end
    PNTSRC=[PNTSRC; ...
        ughb*iPer LRCghb ughb*concPol       ughb*ItypeGhb; ...
        uriv*iPer LRCriv uriv*concRiv(iPer) uriv*ItypeRiv ];
end
%% running scenarios

ICBUND=IBOUND;

save underneath Conf wel cfresh