% nmazex2 (nlagen mazure, example 1)
% PAO 2-3 november 1999, NUMERIEK MODELLEREN, mazure voorbeeld 2
%
% Situatie van Wassen (1990), doorsnede vanaf de Utrechtse Heuvelrug door de Bethunepolder
% tot in de Loosdrechtse plassen
% Theo Olsthoorn, 991014

%Het netwerk
x=[0,5,10,15,17.5,18.5,19,19.5,19.75,20,...
      20.25,20.5,21,22,25,30,40,50,75,logspace(2,4,5)]; dx=diff(x);    Nx=length(x);
z=[0,2,5,10,15,20,25,27,29,30,31,33,35,40,50,60,75,logspace(2,4,8)]'; dz=diff(z);   Nz=length(z);


%Bovenranvoorwaarde, afgelezen uit fig. 2, p84 van Wassen (1990)

k=45*52*ones(Nz-1,Nx-1);
      
% gegeven stijghoogten numerieke model, eerst matrix NaN's dan fixed heads invullen linker en bovenrand
FH=NaN*ones(Nz,Nx); FH(1,:)=0;
FH(10,[1:10])=2.1;

% model runnen, mfile flatmeshctrd staat bij mij in subdirectory "flatmodel"
% tic, toc meet de rekentijd
[Phi,Q]=flatmeshctrd(x,z,k,k,FH,[]);

% plotten van de stijghoogten en de stroomfunctie
figure
contour(x,z',Phi,[0:0.1:2.5]);							% eerst contouren anders zet image plaatje op zijn kop (waarom??)
set(gca,'xlim',[0,50,],'ylim',[0,50]);

Qsleuf=sum(abs(Q(:)))*2
qsleuf=Qsleuf/40


q0=2.1/30*45*52


%Het netwerk
z1=[flipud(-logspace(0,log10(300),5)');z]; Nz1=length(z1);

%Bovenranvoorwaarde, afgelezen uit fig. 2, p84 van Wassen (1990)

k1=45*52*ones(Nz1-1,Nx-1);
      
% gegeven stijghoogten numerieke model, eerst matrix NaN's dan fixed heads invullen linker en bovenrand
FH1=NaN*ones(Nz1,Nx); FH1(1,:)=0;
FH1(15,[1:10])=2.1;

% model runnen, mfile flatmeshctrd staat bij mij in subdirectory "flatmodel"
% tic, toc meet de rekentijd
[Phi1,Q1]=flatmeshctrd(x,z1,k1,k1,FH1,[]);

% plotten van de stijghoogten en de stroomfunctie
figure
contour(x,z1',Phi1,[0:0.1:2.5]);							% eerst contouren anders zet image plaatje op zijn kop (waarom??)
set(gca,'xlim',[0,50,],'ylim',[0,50]);

Qsleuf1=sum(abs(Q1(:)/2))*2
qsleuf1=Qsleuf1/40

fprintf('\n%.0f, %.0f, %.3f, %.3f\n',Qsleuf,qsleuf,Qsleuf1,qsleuf1);