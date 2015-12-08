% GGOR analytic
% Doorsnedegemiddelde grondwaterstand tussen sloten als functie van de tijd
% voor gegeven doorsnede met constant doorlaatvermogen en weerstand rustend
% op een regionale aquifer met vast peil (niet gebruikt) en vaste opwaarts
% positieve gegeven kwel per tijdstap.
%
% Tijdstappen zijn willekeurig maar de randvoorwaarden en inputs moeten
% tijdens elke tijdstap bekend zijn. Daarom wordt standaard gewerkt met 1
% dag. Echter het programma berekent zelf de tijdstaplengte uit de invoer
% reeks van het neerslag overschot.
%
% Wanneer grondwaterstand het maaiveld berekt, vindt oppervlakkige afvoer
% plaats.
%
% Zie ook: numerieke variant gebaseerd op MODFLOW
%
%
% TO 100905 101021 150316
%
% Copyright 2009 2015 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear variables; close all;

basename='HSV';

%% Lees de meteofile
tne = getMeteo('PE-00-08.txt');

%% Lees database in met de eigenschappen van de percelen

data    = dbfread(basename);   % lees database met de percelen

parnams = deblank({data.fieldname}); % verwijder spaties uit veldnamen

% Loop de veldnamen af en pas ze zonodig aan
FID=data( strmatchi( 'FID2', parnams, 'exact')).values;
ZP = data( strmatchi( 'ZP_ACT' , parnams)).values;
WP = data( strmatchi( 'WP_ACT', parnams)).values;
L  = data( strmatchi( 'L', parnams,'exact')).values;
xc = data( strmatchi( 'xc', parnams,'exact')).values;
yc = data( strmatchi( 'yc', parnams,'exact')).values;
k  = data( strmatchi( 'Kh' , parnams, 'exact')).values;
mu = data( strmatchi( 'MU_2', parnams, 'exact')).values;
c  = data( strmatchi( 'CDEK_new' , parnams, 'exact')).values;
q  = data( strmatchi( 'KWELGEMR_2', parnams)).values;
phi= data( strmatchi( 'PHI2', parnams, 'exact' )).values;
AHN= data( strmatchi( 'AHN_MED2', parnams)).values;
AREA=data( strmatchi( 'AREA_METER', parnams)).values;

GLGDBF=data(strmatchi('GLG_FF', parnams)).values;
GHGDBF=data(strmatchi('GHG_FF', parnams)).values;
GVGDBF=data(strmatchi('GVG_FF', parnams)).values;

%% Special actions to make all values compatible and make sense
AREA(AREA<1e3)          = 1e3;% Only use parcels > 1 hectare
q(isnan(q))             = 0.0;
c(c>1000)               = 1000;
L(L<20)                 = 20.0;
b                       = L/2;
AHN( abs( AHN) > 10.0)  = 0.0;
hLR                     =[ZP,WP]; iZP=1; iWP=2;
GLGDBF( GLGDBF < -5)    = NaN;
GHGDBF( GHGDBF < -5)    = NaN;
GVGDBF( GVGDBF < -5)    = NaN;

% Use defaults for missing parameters
D                       = 6.0 * ones(size(b)); % [m] default thickness of top layer
Omega                   = 2.0 * ones(size(b)); % [ m ] default width of ditch
csl_in                  = 1.0 * ones(size(b)); % [ d ] default entry resistance of ditch
csl_out                 = 0.5 * ones(size(b)); % [ d ] default exfiltration resistance of ditch

%% Rounding numbers
rnd{4} = {'phi', 'AHN'};
rnd{3} = {'GP', 'ZP', 'WP', 'mu', 'GLGDBF','GHGDBF','GVGDBF'};
rnd{2} = {'k','c'};
rnd{1} = {'L','xc','yc'};
for i=numel(rnd):-1:1
    for j=1:numel(rnd{i})
        eval([rnd{i}{j} '= round(' rnd{i}{j} ','  num2str(i-1) ');']);
    end
end


%% Check input statistices
parnams = who;
fprintf('\nOverview of database parameters:\n');
for i = 1:numel(parnams)
    eval(['X=', parnams{i} ';']); X=X(:,1);
    if size(X,1)==numel(FID)
        fprintf('%-10s %12.4g %12.4g %12.4g\n',parnams{i},min(X),mean(X),max(X));
    end
end
fprintf('\n');

%% Simulation using analytical solution with constante layer thickness and prescirbed flux

% Add varying ditch level
DV     = datevec(tne(:,1));  month=DV(:,2);
winter                      = true(size(month));
winter(month>=4 & month<=9) = false;
lambda = sqrt( k .* D .* c );
iIntrede= 1;
iUittrede=2;
Sigma(:,iIntrede ) = (b ./ c) .* 2 .* (csl_in  ./ Omega + log(D ./ Omega)/pi ./ k) + b./lambda .* coth(b ./ lambda) -1;
Sigma(:,iUittrede) = (b ./ c) .* 2 .* (csl_out ./ Omega + log(D ./ Omega)/pi ./ k) + b./lambda .* coth(b ./ lambda) -1;
S = ones(size(Sigma(:,1)));

Dt=diff(tne(:,1));
N  = tne(:,2)-tne(:,3);

h  =zeros(numel(L),size(tne,1));  % head in the shallow aquifer
if winter(1)
    h(:,1) = hLR(:,iZP);
else
    h(:,1) = hLR(:,iWP);
end

for it=1:length(Dt);
    
    if winter(it), season=iWP; else season=iZP; end;
    
    I     = ( h(:,it) > hLR(:,season) ); % infiltrating ditches
    S( I) = Sigma( I,1);  % infiltrating ditches
    S(~I) = Sigma(~I,2);  % exfiltrating ditches
    T     = c .* mu .* S;
    E     = exp(-Dt(it))./T;

    h(:,it+1)=hLR(:,season) + ( h(:,it) - hLR(:,season) ) .* E + ...
                 c .* S .* ( N(it) + q ) .* ( 1 - E );

    fprintf('.');
    if rem(it,100)==0, fprintf('\n'); end
    %qditch(:,it)
    %qstor(:,it)
    %qRunoff(:,it)
end
fprintf('\n');

figure; hold on; xlabel('time'); ylabel('NAP [m]'); title('Berekende h');
plot(tne(:,1),h([1:5, end-5:end],:))
datetick()

%% Compute GXG for the analytical simulation

[GLG,GVG,GHG]=getGXG(h(3,:),tne(:,1),0,'plot');

%% Running water balance during the simulation
