%% Aquifer test partial penetration
% Shabbir A.S.Sayed
% Observed Drawdown Pattern Around a Well Partially Penetrating a Vertially
% Extensive Water -Table Aquifer.
% (Groundwater, Vol 22, No. 2, 1984, p 148-153)
%
% PhD Naveed Alam (2012)

%% Problem
% Determine transmissivity from well tests for the Pakistani Punjab
% situation in which the depth of the aquifer is unknown.

%% Approach
% The model is axial symmetric. It contains several rows, each of which is
% a separate model and pumping test. The conductivity in y-direction is set
% to essentially zero to prevent contact between the rows (CHANI in the LAY worksheet).
%
% TO 120924

clear variables;
close all;

%% Parameters for the problem
basename='AquiferTestPP';

AXIAL      = true;
GREP       = 'STRESS';
BACKGROUND = false;

[Dhdr,Data    ]=getExcelData(basename,'data' ,'Vert');
[Whdr,WellData]=getExcelData(basename,'wells','Hor');
[Ghdr,GeoData ]=getExcelData(basename,'geo'  ,'Hor');

%% Data
tstart     = datenum(2013,1,1);
R	       = Data(strmatchi('R'         ,Dhdr,'exact'));

%% Generate all other matrices
ztop       = GeoData(:,strmatchi('ztop',  Ghdr));
zbot       = GeoData(:,strmatchi('zbot',  Ghdr)); 
hk         = GeoData(:,strmatchi('kh',    Ghdr));
vk         = GeoData(:,strmatchi('kv',    Ghdr));
peff       = GeoData(:,strmatchi('peff',  Ghdr));
ss         = GeoData(:,strmatchi('ss',    Ghdr));
sy         = GeoData(:,strmatchi('sy',    Ghdr));

strthd     = GeoData(:,strmatchi('strthd',Ghdr));

%% Grid vertically based on given geology
Z      = [ztop(1); zbot];
xGr    = [0 logspace(log10(0.125),log10(1e5),98)];

%% Generate a grid to grab the geology and the screen length
grOld  = gridObj(xGr,[],Z);          % geology layer extent from 'geo'

%% Arrays for old grid (just one column)
hk     = grOld.const(hk);
vk     = grOld.const(vk);
peff   = grOld.const(peff);
strthd = grOld.const(strthd);
sy     = grOld.const(sy);
ss     = grOld.const(ss);

well = wellObj(basename,'wells');    % well with screen length from 'wells'

%% Vertical grid, refined about the screen tips
ztscr = well(1).z(  1); % screen top
zbscr = well(1).z(end); % screen bot

% vertical grid, refined near top and bottom of screen
zGr    =[sinespace(Z(1),ztscr,15,-pi/2,0), ...
         sinespace(ztscr,zbscr,30,0,pi), ...
         sinespace(zbscr,Z(end),50,0,pi/2)];

%fprintf('%12.2f\n',zGr);

%% Scenario values (we vary only vertical anisotropy and aquifer depth below the screen)
ani  = [1 3 10 30 100];              % vertical anisotropy values
gr   = gridObj(xGr,[],zGr);
ILay = find(gr.zm<well(1).z(end)); % layers below the screen will be varied
ILay = 1; % no depth variants tested

%% Determines the number of options (variants, scnarios) as the product
%  of the number of anisotropy factors and the number of layers below the
%  screen. This translates into the number of rows of the MODFLOW model to
%  allow computing all variants in one model
%  Make sure the horizotnal anisotropy in y-direction in LAY worksheet is zero.

Nrow = numel(ani)*numel(ILay);
yGr= -0.5 + (0:Nrow);

%% Put a well in every row
for iw=Nrow:-1:1
    well(iw)      = well(1);
    well(iw).y    = yGr(iw)+0.5;
    well(iw).nr   = iw;
    well(iw).name = sprintf('well%02d',iw);
end

%% The grid for the entire model, one row per scenarion case
gr     = gridObj(xGr,yGr,zGr,gr.LAYCBD,gr.MINDZ,AXIAL);

%% Transfer geogrid to final grid (one column only)
hk     = gridsTransfer(grOld.zGr,hk    ,gr.zGr,'geometric','z');
vk     = gridsTransfer(grOld.zGr,vk    ,gr.zGr,'geometric','z');
peff   = gridsTransfer(grOld.zGr,peff  ,gr.zGr,'geometric','z');
ss     = gridsTransfer(grOld.zGr,ss    ,gr.zGr,'geometric','z');
sy     = gridsTransfer(grOld.zGr,sy    ,gr.zGr,'geometric','z');
strthd = gridsTransfer(grOld.zGr,strthd,gr.zGr,'geometric','z');

%% Extend the grids into the y direction
HK     = gr.const(hk);
VK     = gr.const(vk);
PEFF   = gr.const(peff);
SS     = gr.const(ss);
SY     = gr.const(sy);
STRTHD = gr.const(strthd);

IBOUND = gr.const(1); % No fixed heads, transient

%% well-bore storage is ignored

%% We will use unconvertible layers for speed and set Ss(:,:,1) accordingly
SS(:,:,1) = SS(:,:,1) + SY(:,:,1)./gr.DZ(:,:,1);

%% Put the wells into the grid, well in every row

well = wellObj(basename,well,gr,HK,'PER');

%% implement scenarios
for iani = 1:length(ani); % for the chosen values of the vertical anisotropy
    % adapt vertical conductivity
    IRow = (iani-1)*length(ILay)+ILay-ILay(1)+1;
    VK(IRow,:,:)=HK(IRow,:,:)/ani(iani);

    for iLay = 1:length(ILay) % for each value of anisotropy
        % inactivate one by one the layers below the screen
        % to simulate incresing aquifer depth
        IBOUND(IRow(iLay),:,ILay(iLay:end))=0;
    end
end

clear grOld
save underneath