% CIE 5440 23 April 2013, Terwisscha asignment
% Adapted to prepare the presentation for NHV 2013-11-27 in Utrecht
% regarding the article by Van den Akker, Stromingen (19)2.

% For this simple case we use a two layer model that is axially symmetric.
% The top layer is the cover layer and the second is the regional aquifer.

% Drains will be used as areal boundaries because they will not allow
% infiltration when the water table sinks to below their elevation.

basename = 'terwisscha';
save name basename

%% Parameters
ayear = 365.24;     % [d/y]

kD1 =   50; D1 =  30; c1 =  100;  % cover layerlayer
kD2 = 4500; D2 = 100; c2 = 2000;  % first aquifer
kD3 = 4500; D3 = 100;            % second aquifer

% scenarios
vdAkker_National = false;

% Defaults
        extDepth = 1.5;  % [m] extension depth in PER sheet
        R  = 3000;   % [m] radius of nature area
        R1 = 8500;
        Nlay=2;
        syNature= 0.25;         % [-] specifid yield in nature area
        syAgri  = 0.10;         % [-] specific yield in aricultural area

[PERnams,PERvals] = getPeriods(basename,'PER');
exDp = mean(PERvals(:,strmatchi('EXDP',PERnams)));
evtr = PERvals(:,strmatchi('EVTR',PERnams));
rech = PERvals(:,strmatchi('RECH',PERnams));
sun  = evtr-rech>0;

scenario = 4;
        
% Special, only what differs from defaults
switch scenario
    case 1 % vdAkker, De Glee, Dupuit Theis
        R  = 1;        
    case 2 % Defaults
    case 3 % Sy in nature and agricltural land
        % make sure extension depth is set to 150 m instead of 1.5 m
        extDepth = 150;
        syNature= 0.10;
        syAgri  = 0.10;
    case 4 % Deep aquifer included
        Nlay = 3;
    case 5 % Model radius at 20000
        R1 = 20000;
    case 6 % Model radius at 2000, deep aquifer included
        Nlay=3;
        R1 = 20000;
end

% Make sure that ET reduction is set correctly for case
if extDepth~=exDp
    error('Change EXDP in PER sheet from %g to %g for scenario %d',exDp,extDepth,scenario);
end


hk      = [kD1/D1  kD2/D2  kD3/D3]; % [m/d]  horizontal conductivit of cover layer

ss      = 1e-5;         % [1/m]; specific elastic storage coeffiient of all layers

Omega   =  1.5;        %  [m] wet circumference of ditches
w       =  0.5;         % [d] entry resistance of ditches

b       =   70;         % [m] half distance between parallel ditches

iWell  =   5;           % [-] zone number for well
iNature=   6;           % [-] zone number nature area
iAgri  =   7;           % [-] zone number agricultural area

% Entry resistance of agricultural land, taking into account both ditches
% and their entry resistance, as is documented in the assignment.
lambda1 = sqrt(kD1*c1);

% Resistane between ditches and phreatic aquifer
cAgri   = 2*b*(w/Omega + log(pi*D1/Omega)/(pi*hk(1))) + c1*((b/lambda1)*cosh(b/lambda1)/sinh(b/lambda1)-1);

%% Grid
Top= 10;

%rGr = [  0 logspace(0,log10(R1),121)];  % grid from 1 to 10000 m
rGr = sinespace(0,R1,121,0,pi/2);
% Top of model is >0, to allow for the free water surface
% the starting heads 0
% So the average thickness of the top layer is still d.
% Make sure LAYCON==1 in sheet LAY to allow for water table behavior
zGr = Top - cumsum([0 D1 D2 1 D3]) ;      % layer elevations

if Nlay<=2
    gr = gridObj(rGr,[-5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5],zGr(1:Nlay+1),'AXIAL',true); % well (y=0) is in first row
else
    gr = gridObj(rGr,[-5.5 -4.5 -3.5 -2.5 -1.5 -0.5 0.5],zGr(1:Nlay+2),'AXIAL',true,'LAYCBD',[0 1]); % well (y=0) is in first row
end

if vdAkker_National
    if Nlay<=2
        gr = gridObj([0 pi*R1^2],gr.yGr,zGr(1:Nlay+1),'LAYCBD',gr.LAYCBD); % not axial
    else
        gr = gridObj([0 pi*R1^2],gr.yGr,zGr(1:Nlay+2),'LAYCBD',[0 1]); % not axial
    end
end

% Boolean to indicate nature and agricultural area
agri       =  gr.Xm>=R;
nature     = ~agri;
drnMdl     = ismember(gr.Ym,[0 -1 -4 -5]);
ghbMdl     = ~drnMdl;
%% Model arrays

IBOUND     = gr.const(ones(Nlay,1)); % no fixed head boundaries
HK         = gr.const(hk(1:Nlay));
VK         = gr.const(hk(1:Nlay)); VK(:,:,1) = D1/2/c1;
VKCB       = gr.constCB(1/c2);
SS         = gr.const(ones(Nlay,1)*ss);
SY         = gr.const(zeros(Nlay,1)); SY(agri) = syAgri; SY(nature) = syNature;
STRTHD     = gr.const(zeros(Nlay,1));

well = wellObj(basename,'wells',gr,HK,{'PER','Q'},'fill',false);

%%
% Boundary conditions for drains
% The conductance CDr of a drain in a cell if the drain is uniformly
% distributed equals CDr = gr.Area/c2;
CDr = NaN(size(gr.AREA));
CDr(agri(  :,:,1)) = gr.AREA(agri(  :,:,1))/cAgri;

zoneValsDRN = {true,  0, CDr(agri & drnMdl)};
zoneValsGHB = {true,  0, CDr(agri & ghbMdl)};

DRN = bcnZone(basename,'DRN',agri & drnMdl,zoneValsDRN);
GHB = bcnZone(basename,'GHB',agri & ghbMdl,zoneValsGHB);

%% De Glee

sDeGlee = well(1).Q(1)/(2*pi*(kD1+kD2))*besselk(0,gr.xm/sqrt(kD2*cAgri));

save underneath b R R1 sDeGlee scenario kD1 kD2 scenario sun