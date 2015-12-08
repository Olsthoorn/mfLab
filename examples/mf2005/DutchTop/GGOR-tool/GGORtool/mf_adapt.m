%% GGOR tool
%
% The GGOR tool models the Dutch geohdyrological top system that consists
% of a shallow Holocene aquifer underlain by a
% semi-pervious layer lying on top of a larger regional aquifer.
% The shallow aquiifer is intersected by many dithces in which the surface
% water level is maintained to manage the groundwater system.
%
% We simulate two parallel ditches in this shallow top system underlain by
% the confining bed and the regional aquifer. The ditches may cut through
% the confining layer.
%
% Input to the model is recharge, seepage from below through the confining
% bed and possibly infiltration from the ditches. Extraction consists of
% evapotranspiration and downward leakage plus discharge to the
% ditches.
%
% The specific yield is constant in this model. Making it dependent is postponed
% to a future extension, if required.
%
% When this water table intersects ground surface, surface runoff occurs.
%
% There is no limit to the infiltration capacity of the soil.
%
% Evapotranspiration may be made depdendent on the extincktion depth.
%
% The exchange between ditch and aquifer passes an entry or exit resistance
% that may be different.
%
% Because we enforce that the managed ditch level is the on both sides of
% the modeled parcel, it suffices to simuulate only half parcel.
%
%
% TO 100905 101021 151112

%clear variables; close all;

basename='GGOR_Aetsveld';  % Polder onder Weesb geplitst door ARK

%% Read meteo time series

tne    = getMeteo('../Meteo/PE-00-08.txt');

if true
    tne(:,3) = 0;

    tt =  365*[0  1  2  3  4  5  6 7];
    pp = 0.01 *[0  1 -1  1 -1  1 -1 -1];

    for it=1:numel(tt)
        tne((180+tt(it)):end,2) = pp(it);
    end

end

%% set wat winter is
DV     = datevec(tne(:,1));        % tne(:,1) is simulation time
winter = DV(:,2)>=10 | DV(:,2)<=3; % winter is a logical vector

%% Maximum parcel width. Parcles > Lmax have subdrainage systems
Lmax = 100; % [m] want percelen > Lmax hebben tenminste greppels die niet op de kaart staan (wel te zien op luchtfoto's)

%% Read the data for all parcels
% dbfRead yield the columns in the dbfFile
% the fieldnames correspond with the properties of the dbfCol

dbfCol  = dbfread(['../GGOR_aetsveld/' basename]);

fldNms   = deblank({dbfCol.fieldname}); %remove blanks from field names

%% Romove all parcels from the database with incomplete dbfCol (having NaNs)
Iparcels = true(numel(dbfCol(1).values),1);

% Mark the dbfColumns with complete numeric information
for ic=1:numel(fldNms)
    if isnumeric(dbfCol(ic).values)
        Iparcels = Iparcels & ~isnan(dbfCol(ic).values);
    end
end

% Only keep the fields with non NaN numeric data
for ic = 1:numel(fldNms)
    dbfCol(ic).values = dbfCol(ic).values(Iparcels);
end

%% Choose a subset of the parcels base for exercising with the GGOR tool

% Skip this if you want to run a full simulaiton using all parcel data
parcels  = 1:50; % set of parcels to consider
Nparcels = numel(parcels);
for ipar = 1:numel(dbfCol)
    dbfCol(ipar).values = dbfCol(ipar).values(parcels);
end


%% Continue with the cleaned up database
% Convert the dbfCol to usable arrays in the Matlab workspace

% Rule to fetch the dbfCol
Fetch  = @(fld) dbfCol(strmatchi(fld,fldNms,'exact')).values;

% Rule to convert a scalar to a vector for all parcels
const  = @(val) val * ones(Nparcels,1);

FID    = Fetch('FID2');

xCtr   = Fetch('XC');     % parcel center x
yCtr   = Fetch('YC');     % parcel center y
ZP     = Fetch('ZP_ACT'); % summer ditch water level 
WP     = Fetch('WP_ACT'); % winter ditch water level
L      = Fetch('L');      % full width of parcel
L(L>Lmax) = Lmax;
L(:)   = 80;              % test to test accuracy
b      = L/2;             % half width of parcel
D1     = Fetch('DDEK');   % thickness of holocene top aquifer (cover layer)
DCB    = const(0.01);     % dummy confining bed (CB) thickness below cover layer (represents its vertical resistance)
D2     = const(30);       % thickness of pleistocene aquifer (2dn layer)
hk1    = Fetch('Kh');     % [m/d] kh in holocene layer
hk2    = const(30);       %[m/d] kh default in regional aquifer
c      = Fetch('CDEK_new');   % [ d ] parcel resistance
vk1    = hk1  ;           % [ m/d ] vk cover layer (in fact dummy, because we use VKCB)
vk2    = hk2/500;         % [ m/d ] default vertical cond regional aquifer
sy1    = Fetch('MU_2');   % parcel Sy
sy2    = const(0.22);     % in case regional aquifer becomes phreatic
ss1    = const(1e-5);     % [1/m] elastic storativity --> optie gebruik dit voor waterberging op maaiveld
ss2    = const(1e-5);     % [1/m] same regional aquifer
q      = Fetch('KWELGEMR_2'); % [m/d] upward seepage
q      = 0 * q;           % TEST
phi    = Fetch('PHI2');       % [ m ] head in regional aquifer (not used)
AHN    = Fetch('AHN_MED2'); AHN(abs(AHN)>10)=0; % ground elevations
AREA   = Fetch('AREA_METER'); % [m2] true parcel area from GIS
wd     = const(2.0);     % [ m ] default width of ditch
drnFac = const(4.0);     %  cIn/cEx;
cEx    = const(5.0);     % [ d ] exfiltration resistance ditch
cIn    = cEx .* drnFac;  % [ d ] default entry resistance of ditch
dd     = const(0.6);     % [ m ] default ditch depth
omega  = wd+2*dd;        % [ m ] default circumference of ditch--> check of dit goed gaat verderop
cdr    = const(1.0);     % [ d ] default drain resistance (as surface resistance)

Nparcel= numel(FID);   % number of parcels

%% Specify ditch levels over time
hDitch = WP * winter' + ZP * ~winter';   % set ditch level

%%  Coordinates for the model grid
dx  = 1;  % cell width choice

% cell boundary coordinates using a tiny left cell to set boundary
% conditions at virtually x = 0. (May be altered layer).
xGr = [ -0.001 0:dx:round(max(b)) ];

% One parcel in each row of the model (set anisotropy to 10-10 in LAY sheet) 
yGr = 0:1:Nparcel;

Nx  = numel(xGr)-1;
Ny  = numel(yGr)-1;

% Rule to convert a parcel vector every model cell in one layer
lay = @(var) var * ones(1,Nx);

%% Z of model grid

deltaZ = 0.001;  % use this for the confining bed at bottom of 1st aquifer

Z   = NaN(Ny,Nx, 3); % allocate memory for Z grid

Z(:,:,1) = lay(AHN);                     % ground surface elevation
Z(:,:,2) = Z(:,:,1) - lay(D1);           % bottom of cover layer
Z(:,:,3) = Z(:,:,1) - lay(D1) - deltaZ;  % bottom of confining bet (=cover layer)
Z(:,:,4) = Z(:,:,1) - lay(D1) - lay(D2); % bottom of regional aquifer

%% Generate grid Object

% guarantee MINDZ not smaller than thickness of (dummy) confining bed
% use LAYCBD to define confining bed at bottom of cover layer
gr = gridObj(xGr,yGr,Z,'MINDZ',deltaZ,'LAYCBD',[1 0]);

%% Id's for set and remember cell properties
iDRN  = 1;    % zone nr of drains in the model
iGHB  = 2;    % zone nr of ditchs in the model
iRIV  = iGHB; % zone nr of river cells in the model (same as iGHB)
iWVP2 = 3;    % zone nr of regional aquifer (where to inject the seepage)
iWEL  = iWVP2;% zone nr where to inject the seepage (same as regional aquifer)

%% Build IBOUND array useing zone numbers to easily identify zones
IBOUND = gr.const(1);
IBOUND(:,2:end,1) = iDRN;    % zone with drains (ground surface)
IBOUND(:,    1,1) = iGHB;    % location of ditch
IBOUND(:,    :,2) = iWVP2;   % second aquifer

% iWEL and iRIV not necessary to set in IBOUND because equal to iWVP2 and iGHB

%% Set active width of each model equal to parcel width b
% This is done by making cells beyond b inactive (IBOUND==0)

for iy=1:gr.Ny
    IBOUND(iy,gr.xm>b(iy),:)=0;
end

% mark which cells are active
active = logical(IBOUND(:,:,1)>0);

%% Build model arrays
HK            = gr.const(NaN);
HK(:,:,1)     = lay(hk1);
HK(:,:,2)     = lay(hk2);

VK            = gr.const(NaN);
VK(:,:,1)     = lay(vk1);
VK(:,:,2)     = lay(vk2);

VKCB          = lay(deltaZ./c);

SY            = gr.const(NaN);
SY(:,:,1)     = lay(sy1);
SY(:,:,2)     = lay(sy2);

SS            = gr.const(NaN);
SS(:,:,1)     = lay(ss1);
SS(:,:,2)     = lay(ss2);

% Notice if laycon == 0 (e.g. to verify GGOR)
% make sure Sy is used, put it in Ss of first layer
[layHdr,layVals] = getExcelData(basename,'LAY','horizontal');
if layVals(1,strmatchi('LAYCON',layHdr))==0
    SS(:,:,1)    = SY(:,:,1)./gr.DZ(:,:,1);
end

%% Initial heads
STRTHD        = gr.const(NaN);
STRTHD(:,:,1) = lay(hDitch(:,1));
STRTHD(:,:,2) = lay(phi);

%% Stress periods to generate RECH and EVTR for MODFLOW
% Use mf2005 because mf2k has a limit of 1000 stress periods

[~,~,NPER]=getPeriods(basename);

% Check that NPER matches length of time series in tne
if size(tne,1) ~= NPER
    error('tne length = %d does not mathc NPER = %d\n',size(tne,1),NPER);
end

%% Recharge MODFLOW specfied here (make sure INRCH == 0 in PER)

% The recharge input is a Ny*Nx*Nper array but can be one of size (1,1,NPER)
% in which case each value is the recharge for the whole model
% for the corresponding stress period.
RECH   = ones(1,1,NPER);
RECH(:)= tne(:,2);

%% Evapotranspiration for MODFLOW specified here
%  EVTR separate from RECH allows reduction of evapotr. in dry spells
EVTR   = ones(1,1,NPER);
EVTR(:)= tne(:,3);

% SURF   = repmat( {AHN - 30.0} , NPER ,1);  % Zie PER sheet
% Make sure the firs inSurf == -2 so that SURF in the PER sheet
% is interpreted as the distance below the top of the model.

%% Drains at ground surface to compute surface runoff
Cdrn   = gr.AREA ./ lay(cdr);   % Drain conductance
Idrn   = find(IBOUND == iDRN);  % Drain locatioins in network
LRCdrn = cellIndices(Idrn,gr.size,'LRC');

%% General Head Boundaries (GHB) in the left-most cell to simulate the ditch
% We need GHB to include entry and exit resistances separately with the
% ditches; we ignore radial resistance for the time being

Cghb   = 0.5*omega./cIn;      % GHB conductance <-- entry resistance
Ighb   = find(IBOUND==iGHB);  % global indices of ditch
LRCghb = cellIndices(Ighb,gr.size,'LRC');

%% RIV to simulate the exit resistance in combination with GHB
Criv   = 0.5*omega.*(1./cEx-1./cIn);
Iriv   = find(IBOUND==iRIV);
LRCriv = cellIndices(Iriv,gr.size,'LRC');

%% Wells to simulate seepage from the bottom aquifer or vice versa
qWEL        = gr.const(0);  % need two layer because well is in layer 2
qWEL(:,:,2) = gr.AREA .* lay(q); % set prescribed seepage in layer 2
Iwel        = find(IBOUND==iWEL);     % Global indices of well
LRCwel      = cellIndices(Iwel,gr.size,'LRC');

%% need a column of ones of given size to multiply scalars with
uDRN = ones(size(Idrn));
uGHB = ones(size(Ighb));
uRIV = ones(size(Iriv));
uWEL = ones(size(Iwel));

%% Specify boundary conditions (stresses) for the MODFLOW model

for iPer = NPER:-1:1
    if iPer == 1
        DRN{iPer} = [iPer*uDRN LRCdrn Z(Idrn) Cdrn(Idrn)];
        WEL{iPer} = [iPer*uWEL LRCwel qWEL(Iwel)];
    else
        DRN{iPer} = [-iPer ones(1,5)]; % as previous
        WEL{iPer} = [-iPer ones(1,4)]; % as previous
    end
        
    if iPer==1 || (winter(iPer) * winter(iPer-1) == false) % change of season
        if winter
            GHB{iPer} = [iPer*uGHB LRCghb WP Cghb];
            RIV{iPer} = [iPer*uRIV LRCriv WP Criv WP];
        else
            GHB{iPer} = [iPer*uGHB LRCghb ZP Cghb];
            RIV{iPer} = [iPer*uRIV LRCriv ZP Criv ZP];
        end
    else
        GHB{iPer} = [-iPer ones(1,5)]; % as previous
        RIV{iPer} = [-iPer ones(1,6)]; % as previous
    end
end

%% Simulation using analytical solution with constant layer thickness and prescibed flux

% Be aware though, that this simulation totally ignores surface runoff
% and evaporation reduction durign dry spells.
% Therefore the numerical and analytical must differ unless the numerical
% settings are set such that the analytical circumstances are guaranteed
% in the numerical model.
% Lastly, the analytical solution works with constant aquifer thickness (to
% avoid zero thickness in case of large downward flow.

kD1     = hk1 .* D1;

Dt      = diff(tne(:,1));
lambda  = sqrt(kD1.*c);
wEx     = D1 .* cEx ./ (0.5 * omega);
wIn     = D1 .* cIn ./ (0.5 * omega);
SigEx   = b./c .* wEx./D1 + b./lambda .* coth(b./lambda) - 1;
SigIn   = b./c .* wIn./D1 + b./lambda .* coth(b./lambda) - 1;

Tex     = sy1 .* c .* SigEx;  % Time constant when exfiltrating
Tin     = sy1 .* c .* SigIn;  % Time constant when infiltrating
cSigEx  = c .* SigEx;         % factor when exfilrating
cSigIn  = c .* SigIn;         % factor when infiltrating

%% Analytical simulation

% Allocate space to store the analyical heads (average head in cross sections)
hAnalytic      = NaN(Nparcel,NPER);  % head in the shallow aquifer
hAnalytic(:,1) = hDitch(:,1);        % initialize first values

% Simulate while diffentiating when in- and exfiltration occurs
% All sections are simulated simultaneously
for it=1:length(Dt);
    Lex = hAnalytic(:,it)>hDitch(:,it);
    hAnalytic(:,it+1) = hDitch(:,it) + ...
        Lex .* ((hAnalytic(:,it)-hDitch(:,it)).*exp(-Dt(it)./Tex)...
             + cSigEx.*(RECH(it)-EVTR(it)+q).*(1-exp(-Dt(it)./Tex))) + ...
       ~Lex .* ((hAnalytic(:,it)-hDitch(:,it)).*exp(-Dt(it)./Tin)...
             + cSigIn.*(RECH(it)-EVTR(it)+q).*(1-exp(-Dt(it)./Tin)));

    % Let us hear from you ...
    fprintf('.'); if rem(it,100)==0, fprintf('\n'); end
end
fprintf('\n');

%% Analytically compute the final head for steady state
% This is for checking GGOR, may be skipped here, see steadySolution.m

Nend = RECH(end) - EVTR(end); % Recharge in last time step

% Ditch resistance in last time step
w = wEx .* (hAnalytic(:,end)>hDitch(:,end)) + wIn .* (hAnalytic(:,end)<=hDitch(:,end));

% x-coordinate along the analytical cross section
x  = bsxfun(@minus, b, gr.xm);

% Analytical steady state solution [m]
hx = bsxfun(@times, phi + Nend.*c, active) - ...
     bsxfun(@times, (phi +Nend.*c-hDitch(:,end)) ./ ...
         ((lambda./c .* w./D1).*sinh(b./lambda) + cosh(b./lambda)), ...
         cosh( bsxfun(@times, x, 1./lambda) ) );

% Analytical steady state solution for flow to or from the ditch [m2/d]
Q = (phi + Nend.*c - hDitch(:,end)) .* sinh(b./lambda) ./ ...
    ((w./D1).*sinh(b./lambda) + (lambda./kD1) .* cosh(b./lambda));

%% Save whatever non-modflow data is necessary in mf_analyze
save underneath FID tne phi q kD1 D1 w c b omega lambda hDitch hAnalytic AREA
