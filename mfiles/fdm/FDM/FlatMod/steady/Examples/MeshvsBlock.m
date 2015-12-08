%PATO cursus 2 en 3 november 1999, Numeriek modelleren
%meshvsblock (nlagen mazure, example 1)
%Vergelijking van het vlakke mesh-centred model met het vlakke block-centred model
% T.N.Olsthoorn 991029

clc; for i=1:30, fprintf('\n'); end			% clear screen, cursor down 30 lines

kD=[1000;500; 750];
c =[ 400;800;1100];
Phi0=[-2;-1;-3];

% numerieke oplossing n-lagenprobleem:
dz=[5;20;20;35;10;50]; z=[0;-cumsum(dz)];
x=[logspace(0,log10(100001),41)]-1;
x=[0:50:3000];

%laagdikten en zwaarden van de knooppunten
k =[dz(1)/c(1);kD(1)/dz(2);dz(3)/c(2);kD(2)/dz(4);dz(5)/c(3);kD(3)/dz(6)];
kx=k*ones(1,length(x)-1);
kz=kx;
FH=NaN*ones(length(z),length(x)); FH([2,4,6],1)=Phi0; FH([3,5,7],1)=Phi0; FH(1,:)=0;
[PhiMC,QMC]=flatmeshctrd(x,z,kx,kz,FH,[]);


% blockcentred
kz(1,:)=0.5*kz(1,:);
FH=NaN*ones(size(kx)); FH([2,4,6],1)=Phi0; FH(1,:)=0;
[PhiBC,QBC]=flatblockctrd(x,z,kx,kz,FH,[]);

fprintf('Contouren van stijghoogten berekend met resp. het mesh- en het block-centred model...\n');
fprintf('Eerst het mesh-centred model...\n');
fprintf('Any key to continue...\n'); pause

figure
patch([x(1),x(end),x(end),x(1)],[z(1),z(1),z(2),z(2)],'m'); hold on
patch([x(1),x(end),x(end),x(1)],[z(3),z(3),z(4),z(4)],'m');
patch([x(1),x(end),x(end),x(1)],[z(5),z(5),z(6),z(6)],'m');

contour(x,z,PhiMC,[-3:0.25:0]);
xlabel('x [m]'); ylabel('z [m]'); title('Doorsnede met 30 stijghoogte- en 30 stroomlijnen');

fprintf('Fraaie isohypsen, ondanks dat we slechts 6 modellagen hebben gebruikt.\n');
fprintf('Any key to continue...\n'); pause

fprintf('Toevoegen van de isohypsen van het blok-centred model..\n');
fprintf('Any key to continue...\n'); pause

xM=0.5*(x(1:end-1)+x(2:end))-0.5*(x(2)-x(1));
zM=0.5*(z(1:end-1)+z(2:end));
contour(xM,zM,PhiBC,[-3:0.25:0]);

xlabel('x [m]'); ylabel('z [m]'); title('Doorsnede met 30 stijghoogte- en 30 stroomlijnen');

fprintf('De isohypsen van het blok-centred model stroken niet met het werkelijke verloop\n');
fprintf('binnen elke afzonderlijke laag. Dit is altijd zo bij deze modellen en dus ook\n');
fprintf('in (VISUAL) MODFLOW.\n');
fprintf('De resultaten zijn echter aan elkaar gelijk omdat de isohypsen van de modellen\n');
fprintf('De resultaten zijn aan elkaar gelijk omdat de isohypsen van de modellen\n');
fprintf('elkaar exact in het midden van elke laag snijden!!\n');
fprintf('Any key to continue..\n'); pause

fprintf('Vergelijken van de berekende stijghoogte met het mesh-centred en block-centred model..\n');
fprintf('Eerst verloop berekend met het mesh-centred model...\n');
fprintf('Any key to continue..\n'); pause

figure;
xLim=[0 1000];
a(1)=subplot(3,1,1); plot(x,PhiMC(2,:),'r'); set(gca,'xLim',xLim); ylabel('head aqf 1 [m]');
title('Stijghoogten in de drie aquifers berekend met mesh- en blok-gecentreerd model');
a(2)=subplot(3,1,2); plot(x,PhiMC(4,:),'r'); set(gca,'xLim',xLim); ylabel('head aqf 2 [m]');
a(3)=subplot(3,1,3); plot(x,PhiMC(6,:),'r'); set(gca,'xLim',xLim); ylabel('head aqf 3 [m]'); xlabel('x [m]');

fprintf('Toevoegen de berekende block-centred model..\n');
fprintf('Any key to continue..\n'); pause

line(xM,PhiBC(2,:),'parent',a(1),'color','b','linestyle','none','marker','+');
line(xM,PhiBC(4,:),'parent',a(2),'color','b','linestyle','none','marker','+');
line(xM,PhiBC(6,:),'parent',a(3),'color','b','linestyle','none','marker','+');

fprintf('U ziet er is geenverschil tussen de uitkomsten van de twee modellen!!\n');