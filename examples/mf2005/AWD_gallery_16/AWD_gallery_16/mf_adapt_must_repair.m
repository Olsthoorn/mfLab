%% 3D MODEL OF Drain 16, to simulate its clogging and cleaning
%
%% Background
%
% mfLab model generation file for MODFLOW
% addtional files needed (see load commands in this file and getExcelData)
%
% DRAIN.mat (shapes of the drain in the Amsterdam Water Supply Dunes
% CANAL.mat (shapes of the canals in the Amsterdam Water Supply Dunes
% POND.mat  (shapes of the recharge ponds)
%
% Drain16.xls    % parameter fie for MODFLOW 
% Dr16_input.xls % additional data (piezometer and waterbodies)
%
% mf_adapt   (m-file (script) to generate the model arrays
% mf_analyze (m-file (script) to visualize the results
% the mfLab environment see
%  http://code.google.com/p/mflab
%

%%
% TO 100210-110813 Recharge basins and Canals surrounding drain 16 and drain 16 itself
% 3D model for the area around Drain 16 in the Amsterdam Water Supply Dunes
% The model includes surrounding recharge basins, extraction canals and a
% seepage pond that fluctuates along with the groundwater but shortcuts
% water horizontally.
% THe canal with the low level in the north (Huppel Canal) attracts water
% from the adjacent recharge basins, and could potentially cause short
% circuit through the drain. However, the model as such does not seem to
% yield this shortcut easily. The shortcut is even more unlikely if the toe
% of the north drain is clogged, as it is likely the case if short cuts
% have taken place through the drain in the past.
% Simulaton with cloggin include is a non-trivial task because of the lack
% of sufficient data of the clogging resistance inside the drain and in its
% gravel pack and slots alongside the drain.

%% What to do next?
%
% Check the model further
% --- dry cells for those where Phi is below top of cel (done but not
%  checked properly)
% --- excact conceptual approach for implemetnation of recharge pond bottom
% resisistance. Same for canals and seepage pond. What to do with the cells
% on top of the canals and ponds? What to do wit their rims. How to use
% inactive cells?
% If ok. Then compare with measured heads inside and outside the drain.
% Make sure that the flow matches the measured seepage. This would allow
% setting the outside and inside resistances through calibration. Wonder if
% a unique result can be obtained with the available data and its
% uncertainty.
% No pipe friction has been included yet.
%%
% TO 100215

%% Notice
% The model uses data DR16.mat, and modelinput.xls
% use these data and make sure it's sufficient and accurate don't hardwire
% inside the code !!

%% Gnerating the input
clear variables; close all

basename = 'Drain16';
verbose = 1; % more output along the way

%AFTERMFSETUP='mf_analyze';

%% Overall model parameters 

kh     = 12; kv=kh;      % [m/d] conductivities
sy     = 0.24;           % [ - ] specific yield
ss     = 1e-5;           % [1/m] specific elastic storage coefficient
havg   = 7.0;            % [ m ] average head for startheads        
KSURFWATER = 1000;       % [m/d] surface water conductivity used in seepage pond
rw     = [0.3   0.3];    % [ m ] drain radius
Rw     = [0.4   0.4];    % [ m ] drain outer radius (outside gravel pack)
c0     = [0.02 0.02];    % [ d ] drain resistance at heel and toe resp
cL     = [0.4   0.4];    % [ d ] drain resistance at heel and toe resp
lambda = [300 300];      % [ m ] drain resistance shape factor
dphi   = [0.4 0];        % [ m ] phi(toe)-phi(heel) drain 1
kDrLong= [100 100];      % [m/d] k of drain in longitudinal direction, this should match the actual
zheel  = 0;              % [ m ] heel or drain elevation
z0     = 0;              % [ m ] drain elevation, all drains in the dunes are about zero msl (NAP)

%% time with observations, choose 1 to set boundary conditions
times(1)=datenum(2009,11,27,12,0,0); % 11/27/09
times(2)=datenum(2009,12,05,12,0,0); % 12/05/09
times(3)=datenum(2009,12,09,12,0,0); % 12/09/09
times(4)=datenum(2010,06,07,12,0,0); % 06/07/10

time=times(4);  % time for this (steady state siimulation

%% Get piezometer basic data and put in a struct
[~,piezprops,txthdr,pieznams]=getExcelData('Dr16_input','Piezometers','Hor');

Piez=piezometerObj(pieznams,piezprops);

%% Add measured heads to Piez

[~,data,~,pieznams]=getExcelData('DR16_input','Heads','H');

for i=1:length(Piez)
    I=strmatchi(Piez(i).name,pieznams,'exact');
    if I(1)>0
        Piez(i).t=excel2datenum(data(I,1));
        Piez(i).h=data(I,2);
    end
end

%% Drains (in our case the of the U11 outlet)
%  Whatever the coordinates in the database, we take the outmost two
%  coordinate pairs and draw a staight line between them, which is the line
%  of the two drains


load('Drain');   % All drain shapes of the dune area

I=strmatchi('Drain_U11',{DRAIN.NAAM});  % only need U11

% The shaft location is one of the piezometer locations. It is the center
% of the two drains, the heels of both.
i=strmatchi('shaftM',{Piez.name});      % get corresponding piezometer

% SHAFT           LL point of drain       UR point of drain
xshaft=Piez(i).x; xm = min([DRAIN(I).x]); xM = max([DRAIN(I).x]);
yshaft=Piez(i).y; ym = min([DRAIN(I).y]); yM = max([DRAIN(I).y]);

%plot(xshaft,yshaft,'yp');
% move shaft exactly onto the line between the two outer most drain points
[xshaft,yshaft]=point2line([xm xM],[ym yM],xshaft,yshaft);

%%
% Generate 2 drain objects (north and south for outer points to shaft = heel
Drain(2)=drainObj('DR16N',[xshaft,     xM], [yshaft,     yM],[z0 z0], [], xshaft, yshaft, zheel);
Drain(1)=drainObj('DR16S',[    xm, xshaft], [    ym, yshaft],[z0 z0], [], xshaft, yshaft, zheel);

%% Lets move the piezometes above the drain exactly above it 

% figure; hold on;
% plot(Drain(1).x,Drain(1).y,'b','linewidth',2);
% plot(Drain(2).x,Drain(2).y,'m','linewidth',2);
% plot([Piez.x],[Piez.y],'ro');

CP=105:116;  % numbers of the C-wells connected with DR16N and DR16S

for i=1:numel(CP) % For all C-wells of DR16N and DR16S
    cwell=sprintf('C%d',CP(i));             % name of C-well (string)
    j=strmatchi(cwell,{Piez.name},'exact'); % look it up in piezometer list
    if j(1)>0  % if found
        % see how much it needs to be move to be exactly above DR16
        [~,~,dx,dy]=point2line([xm xM],[ym yM],Piez(j).x,Piez(j).y);
        
        % next look for related wells C105M_A, C105M_B C105M_C et
        J=strmatchi(cwell,{Piez.name});
        for k=J  % all of them will be move as well
%             plot(Piez(k).x,Piez(k).y,'k');
%             plot([Piez(k).x,Piez(k).x+dx],[Piez(k).y,Piez(k).y+dy],'r-');
%             text(Piez(k).x,Piez(k).y,Piez(k).name);
%                  
            % over the same distance as the C-well it self. So that the
            % wells closest to the drain are at their accurate relative
            % locations with respect to the drain itself.
            Piez(k)=Piez(k).shift(dx,dy);
%             plot(Piez(k).x,Piez(k).y,'bo');
            
        end
    end
end

%% With all piezometers in place we get the waterbodies

load('Canal');NC =numel(CANAL); % Canal data (shapes) from disk
load('Pond'); NP =numel(POND);  % Recharge pond data (shapes) from disk
              NSP=1;            % We also have one seepage pond

% Generate a waterbody object for all of them (memory preallocation)             
waterbody(NC+NP+NSP)=waterBodyObj();

waterbody( 1:NC )       = waterBodyObj(CANAL,'canal'); % Get the canals
waterbody(   NC+(1:NP)) = waterBodyObj(POND, 'pond');  % Add the Ponds

% Seepage pond becomes the last waterbody. We wil generate it seperately
% using the shape digitized in Google Earth. We get the coordinates form
% the kml path file directly, and tranform them into the Dutch coordinate
% system using wgs2rd(EN] directly and adding this to the waterbody.
waterbody(end)=waterbody(end).setKML('DR16SpgPond.kml');
waterbody(end).type    = 'SpgPond';
waterbody(end).gauge   = 'SpgPond';
waterbody(end).name    = 'Seepage Pond Near C105 Drain 16';

%% Connect each water body with its gauge, which are also in the piezometer
% data set. We check the water bodies one by one and look which piezometer
% lies within the waterbody contour:

figure; hold on

% I remembers the water bodies with a piezometer in their contour. The
% other waterbodies are too far away to be relevant for out local model
I=zeros(size(waterbody));

for i=1:numel(waterbody)
    for j=1:numel(Piez)
        if inpolygon(Piez(j).x,Piez(j).y,waterbody(i).x,waterbody(i).y)
            % if its in, we plot the waterbody and its piezometer
            plot(waterbody(i).x,waterbody(i).y,'b');                     % debug
            plot(Piez(j).x,Piez(j).y,'ro');  % debug
            
            I(i)=1; % and remember this waterbody
            % while connecting it to its gauge (in the piezometer set)
            waterbody(i).gauge=Piez(j).name;
            
            fprintf('waterbody %d matches piezometer %d\n',i,j)
        end
    end
end

% finally, do away with all irrelevant waterbodies
waterbody=waterbody(I==1);

%% Get data pertaining to named water bodies from accompanying Excel workbook
% DR16_input.xls. These data are in the worksheet Waterbodies. The data
% is the elevation and hydraulic resistance of the bottoms. We don't use
% x_ctr and y_ctr, which we compute ourselves from the waterbody
% circumference.

[hdr,Data,txtdhr,txt]=getExcelData('DR16_input','WaterBodies','H');

% Find the columns containing the data in the table
iC_entry = strmatchi('c_entry',hdr);
iz_bot   = strmatchi('z_bot',  hdr);
iDescr   = strmatchi('Descrip',txtdhr);

% First deal with the two drains
for i=1:numel(Drain)
    j=strmatchi(Drain(i).gauge,txt(:,1));
    if j(1)>0 && ~isempty(Drain(i).gauge)
        Drain(i).c_entry=Data(j,iC_entry);
        Drain(i).z_bot  =Data(j,iz_bot);
        Drain(i).name   =txt {j,iDescr};
    end
    
    % Add other data, stored in the head of this m-file
    k=strmatchi(Drain(i).gauge,{Piez.name});
    Drain(i).c0      = c0(i);     % resistance of drain heel
    Drain(i).cL      = cL(i);     % resistance of dain  toe 
    Drain(i).lambda  = lambda(i); % resistance shape factor
    Drain(i).dphi    = dphi(i);   % Phi(toe)-Phi(heel)
    Drain(i).rw      = rw(i);     % drain inner radius
    Drain(i).Rw      = Rw(i);     % drain gravel pack raidus
    Drain(i).kDrLong = kDrLong(i);% drain longitudinal conductivity
end

%% Same thing for all water bodie, canals and ponds
for i=1:numel(waterbody)
    j=strmatchi(waterbody(i).gauge,txt(:,1));
    waterbody(i).havg = havg;
    if j(1)>0
        waterbody(i).c_entry=Data(j,iC_entry);
        waterbody(i).z_bot  =Data(j,iz_bot);
        waterbody(i).name   =txt {j,iDescr};
    end
end

%% Rotate and shift the world
% Next we align the world along the model grid and use the central
% shaft as the center of the model universe and roates all objects in the
% model around the central shaft such that the drain aligns with the
% x-axis.

% First get the rotation angle of the drain with respect to the x-axis
alfa=-atan2(diff(Drain(1).y),diff(Drain(1).x)) * 180/pi;

%% Rotate the drain, Canal, Poinds and Piezometers around the shaft.
for i=1:length(    Drain),     Drain(i) =     Drain(i).rotate(xshaft,yshaft,alfa); end
for i=1:length(waterbody), waterbody(i) = waterbody(i).rotate(xshaft,yshaft,alfa); end
for i=1:length(     Piez),      Piez(i) = Piez(     i).rotate(xshaft,yshaft,alfa); end

%% Plot new model in rotated and centered coordinates, to verify if all went well

figure; hold on

% See how nicely each object knows how to plot itself (using dot notation)
for i=1:length(Drain),      Drain(    i).plot('b-',2); end
for i=1:length(Piez),       Piez(     i).plot('ro');  end
for i=1:length(waterbody),  waterbody(i).plot('r-',1,'canal'); end
for i=1:length(waterbody),  waterbody(i).plot('b-',1,'pond');  end


%% The model grid

% Generate model grid coordinates, making sure we can have detail where we
% need it, i.e at steep gradients (ends of the drain and vicinity of drain.

% Choose a suitable bounding box
Bx=[-450 -300 450 300]; W=diff(Bx([1 3])); V=diff(Bx([2 4])); % Model boundary box

% Tell how detailed you want the model to be
Nx1=20; Nx2=60; Nx3=15;       % used in grid generation below
Ny= 60;                       % used in grid generation below

% Select a list of layer elevations thought useful and suitable for this model
zGr=[7 5 3 2.25 1.5 1.0 0.5 0.4 -0.4 -0.5 -1 -1.5 -2.25 -3 -4 -7 -11.9 -12];

% Generate smooth y-grid coordinates using mfLab's sinespace function. This
% allows a very flexible mesh generation. yGr will run from the bottom  to the
% top of the bounding box usign Ny/2 gid lines in both halvs
% end refining near the tips of the two drains. Our smallest cell will
% runform -0.5 to +0.5 m.
yGr=[sinespace(Bx(2) ,-0.5  ,Ny/2,pi/2,0) -0.4 0.4 sinespace(0.5,Bx(4),Ny/2,0,pi/2)]; 

% Generate smooth x-grid coordinates, using drain coordinates. Start left
% at the bounding box to the toe of the left drian, usign sinespace(). Then
% we continue along the drain using sinespace to obtain a fine spaceing
% near the drain tips. Same thing for the right half of the model.
xGr=[sinespace(Bx(   1)       ,Drain(1).x(  1),Nx1,pi/2,0), ...
     sinespace(Drain(1).x(  1),Drain(2).x(end),Nx2,0,pi),   ...
     sinespace(Drain(2).x(end),Bx(     3)     ,Nx3,0,pi/2)];

 % ensure there is a grid line exactly through the shaft. So that we can
 % truly divide the model in a left and a right half the drain with equal
 % size.
xGr(abs(xGr)<2)=0;

% Finally generate a grid object to store any grid information we need.
grid=gridObj(xGr,yGr,zGr);

% show the grid:
grid.plot('xy','c');

%% Add Waterbodies and zones to the model

% We still have to merge the waterbodies with our grid. Here we do this
% with all water bodies, one oafter the other.

I=ones(size(waterbody));  % Remember those who intersect the grid
for izone=1:numel(waterbody)
    % Here we merge, using the grid as argument to the function merge
    % which is a memthod of the waterbody objects. It adds grid info
    % to the waterbody objects, so it then knows where it is in the grid.
    waterbody(izone)=waterbody(izone).merge(grid,izone); % merge
    if isempty(waterbody(izone).I), I(izone)=0; end      % remember
end 

% Do away with the waterbodies that do not intesect the grid.
waterbody=waterbody(I==1);

%% Add the drains to the model

for i=1:numel(Drain)
    Drain(i)=Drain(i).merge(grid); % merge drain into model
end

%% MODLFOw models arrays using grid to set their size

IBOUND= gr.const(1);;  % Boundary condition indicator array

HK     = gr.const(kh);  % Horizontal hydraulic conductivity
HKY    = HK;                  % ky
VK     = gr.const(kh);  % same
SY     = gr.const(sy);
SS     = gr.const(ss);
STRTHD = gr.const(havg); % fixed heads initiation

%% possibly contour the zones to show what we've got

% So let us put the location of the waterbodies as zones into the grid
for i=1:numel(waterbody())
    IBOUND(waterbody(i).I)=waterbody(i).izone;
end

if verbose~=0  % may be skipped by setting verbose to 0
    for iz=1:grid.Nz,
        % if layer not uniform (constant) then contour and plot it
        if ~all(IBOUND(:,:,iz)==IBOUND(1,1,iz))
            contourf(grid.xm,grid.ym,IBOUND(:,:,iz));
        end
    end
end

%% Fix the head in all waterbodies Set the fixed heads of all waterbodies
% We may just as wlel use GHB RIV and DRN to this end, it's a choice.

for i=1:numel(waterbody),
    % Set HK to very in the water phase of the waterbodies. This will
    % ensure almost constant head whithin canals. Altenratively use the
    % fixe heads or the lake package for more sophistication.
    HK(    waterbody(i).I) = KSURFWATER;
    
    % Don't fix the head of all waterbodies except that of the seepage pond 
    if ~strcmpi(waterbody(i).gauge,'spgpond')
        % This si done by making the index in IBOUND negative.
        IBOUND(waterbody(i).I) = -waterbody(i).izone;
    end
end

%Because we have a different KH compared to KV, we apply anisotropy by
%specifying its horizontal valus throught the anisotropy array.
HANI=HKY./HK;
   
%% Other boundary conditions

% Recharge is not so relevant for this study, so we do not include it for now
% Same for evapotranspiration.

% Get stress period information
[pernams,pervals]=getPeriods(basename);
NPER=size(pervals,1);

% we may use different boundary conditions simultaneously. The advantage is
% that the outcomes are easily recognized in the buget file. However we are
% free to use them or not.

RIV=[];  % we use the river boundary for the head in the canals
GHB=[];  % we use the general head boundary for the head in the recharge pomds
CHD=[];  % we use this to set time varying head boundary conditions.
WEL=[];  % we use this for wells (WEL package) if there are wells.

% Generate the lists from which we write the MODLFOW input files
for i=1:length(waterbody)
    
    % Set vertical hydrauci VCONT (=1/c) of the bottom of each of the
    % recharge ponds. See c_entry and their location (given their
    % circumfernce and thei bottom depth z_bot. This VCONT is used by all
    % general head type of boundaries.
    waterbody(i)=waterbody(i).setCond(grid);
    u=ones(size(waterbody(i).I));
    for iPer=1:NPER
        switch waterbody(i).type
            
            % Boundary condition lies for the canals
            case 'canal'  % we use GHB for canals
                GHB=[GHB;...
                    iPer*u waterbody(i).LRC u*waterbody(i).head(Piez,time) waterbody(i).Cond ...
                    ];
                
            % Boundary condition list for the recharge ponds
            case 'pond'  % we use RIV for ponds
                RIV=[RIV;...
                    iPer*u waterbody(i).LRC u*waterbody(i).head(Piez,time) waterbody(i).Cond, u*waterbody(i).z_bot ...
                    ];
        end
    end
end

%% Horizontal flow barrier around every waterbody

HFB=[];
for i=1:length(waterbody)
    HFB=[HFB; ...
       setHFB(grid.xm,grid.ym,1:waterbody(i).izbot,waterbody(i).x,waterbody(i).y,waterbody(i).c_entry)];
end

%% Inside the drains we use CHD as boundary. These are fixed heads which may vary
%  within the stress periods. Note that the well entry ressitance is
%  transferred to the cells in which the drain resides. Therefore, this
%  head is to be considered the head inside the drain. But no actual flow
%  is simulated inside the drain yet. However, this can be included using
%  the MNW (multinode well) package the stream routing package and the CFP
%  (conduit flow package). The lattes is most sophisticated and more
%  compliex.

ip=strmatchi('C110',{Piez.name},'exact');

for i=1:length(Drain)
    u=ones(size(Drain(i).I(:)));
    for iPer=1:NPER
        CHD=[CHD;...
             iPer*u Drain(i).LRC...
             Drain(i).head(Piez(ip),time)'...
             Drain(i).head(Piez(ip),time)' ...
            ];
    end
end

% We unpack the grid a bit, because mf_setup needs it to generate the input
% files for MODFLOW
Z=grid.Z; Dx=grid.Dx; Dy=grid.Dy; Dz=grid.Dz;

%Finally save some extra information that is usefull for postprocessing.
save underneath Piez waterbody Drain grid

%The input for MODFLOW will be generated and MODLFOW will berun when you
%type 
%
%mf_setup<<RETURN>>
%
% After MODFLOW has finishwed type
%
%mf_analyze
%
% To see the results.