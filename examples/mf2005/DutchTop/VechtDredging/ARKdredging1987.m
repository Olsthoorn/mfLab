% ARK canal dredging: Check report by Steenhuis (1987)
% INVLOED van de ZANDWINNING op de KWEL uit het AMSTERDAM-RIJNKANAAL
% (Eng: Influence of sand dredging on the seepage from the Amsterdam-Rhine Canal)
% Report:          DBW/RIZA nota 97.048
% Author:          B.P.C. Steenkamp, oktober 1987
% Publisher:       Dienst Binnenwateren/RIZW, Netherlands
% Main department: Hoofdafdeling Watersystemen
% Departments:     Afdelingen Delta en Meren (WS)

% TO 001203 JB+TO 101123

% Check using obs well/piezometer 31E0166 screen 1 NAP -28.9-29.9 and screen 2 NAP -8.9 -9.9 m


clear variables;
close all

%% Data from report:
XP2l=-5000; XP1l=-72;XSl=-70;XRl=-50;XBl=-30;XBr=30;XRr=50;XSr=70;XP1r=72;XP2r=2450;
HP2l=-6.6;HP1l=-1.9;HSl=-1.9;HDl=-0.4;HRl=-0.4;HB=-0.4;HRr=-0.4;HDr=-0.4;HSr=-1.7;HP1r=-1.7;HP2r=-1.2;
CP2l=100;CP1l=250;CSl=250;CDl=1e6;CRl=250;CB=250;CRr=250;CDr=1e6;CSr=250;CP1r=250;CP2r=250;
kd=4396; vAniso=3;

%% Choose a scenario i.e. a year
scenario=-1981;

switch scenario
    case 1980,  NLay=1; cARK=[250 250 250]; hPeilbuis=-1.75; xPeilbuis=-70;  % voor het baggeren in 1981, meot zijn NAP -1.75 m
    case 1981,  NLay=1; cARK=[ 40 6.5  40]; hPeilbuis=-1.05; xPeilbuis=-70;  % direct na het baggeren   , moet zijn NAP -1.05 m
    case 1987,  NLay=1; cARK=[ 50 40   50]; hPeilbuis=-1.50; xPeilbuis=-70;  % direct na het baggeren   , moet zijn NAP -1.05 m
    case 1993,  NLay=1; cARK=[ 70  70  70]; hPeilbuis=-1.62; xPeilbuis=-70;  % tien jaar na het baggeren, moet zijn NAP -1.62 m
    case -1980, NLay=5; cARK=[250 250 250]; hPeilbuis=-1.75; xPeilbuis=-70;  % voor het baggeren in 1981, meot zijn NAP -1.75 m
    case -1981, NLay=5; cARK=[ 40 6.5  40]; hPeilbuis=-1.05; xPeilbuis=-70;  % direct na het baggeren   , moet zijn NAP -1.05 m
end


Z=[-2 -8 -50]';  Daquitard=abs(diff(Z(1:2))); Daquif=abs(diff(Z(2:3)));

d=1;
D1=ones(NLay,1)*d; D2=ones(NLay,1)*Daquif/NLay; D2(2:end)=D2(2:end)-d;
D=[D1 D2]'; D=D(:); D(1)=Daquitard;

z=[Z(1); Z(1)-cumsum(D)];

x=    [   XP2l   XP1l   XSl   XRl   XBl  XBr   XRr   XSr   XP1r   XP2r    ];	   % zone boundaries
h=    [HP2l   HP1l   HSl   HDl   HRl   HB   HRr   HDr   HSr    HP1r   HP2r];	   % prescribed heads
n=0 * h;				                                                           % precipitation-evapotranspiration
C=    [CP2l   CP1l   CSl   CDl  cARK(1) cARK(2) cARK(3) CDr   CSr   CP1r    CP2r]; % hydraulic resistance of topsystem

kD=kd/Daquif*D2*ones(size(h));

c= ((D1+D2)*ones(size(h))).^2./kD; c(1,:)=C;  c(c<1)=500;

QZ=   zeros (NLay,length(h)-1);

X=-7000:1:7000;
xp=x(x>X(1) & x<X(end));
X=unique([X,xp-0.1,xp+0.1]);

%% print section data

[~,phi,qS]=nsecn(x,kD,c,h,n,QZ,x);  % analytical multi-layer solution

diff(qS)

fprintf('%10s\t%10s\t%10s\n','x','Qx','stijgh');
fprintf('%10s\t%10s\t%10s\n','m hart ARK','m2/d','m +NAP');

for i=1:length(x)
    fprintf('%10.0f\t%10.3f\t%10.2f\n',x(i),qS(i),phi(i));
end

%% analytical multi-layer solution:
[Section,phi,q,s,X,x,kD,c,H,n,Q]=nsecn(x,kD,c,h,n,QZ,X);

shownsecn(phi,q,s,z,h,X,x); % plot result of this solution

ylim=get(gca,'ylim'); set(gca,'ylim',[ylim(1),0]);
plot(X,phi(1,:),'r');

%% numerical solution using finite difference method

% Grid
x=x(x>-Inf & x<Inf);

xGr=X;
zGr=z;

[xGr,yGr,xm,ym,Dx,Dy,Nx,Ny]=modelsize(xGr,z);

FH=NaN(  Ny,Nx);   % prescribed heads
FH(1,:)=       h(1);

FQ=zeros(Ny,Nx);   % prescribed flows

kx=zeros(size(FH)); % horizontal conductivities
kx(1:2:end,:)=Dy(1:2:end)./c(:,1)   * ones(size(xm)); kx(1,:)=kx(1,:)/2;
kx(2:2:end,:)=kD(:,1)./Dy(2:2:end)  * ones(size(xm));

% adapt
for i=1:length(x);
  Ix=find(xm>x(i));
  FH(1      ,Ix)=      h(i+1);  % match with prescribed zones heads
  kx(1:2:end,Ix)=Dy(1:2:end)./c(:,i+1)/2      * ones(size(Ix));
  kx(2:2:end,Ix)=      kD(:,i+1)./Dy(2:2:end) * ones(size(Ix));
end 

%figure; plot(xm,FH(1,:));
%figure; plot(xm,kx(1,:));

IBOUND = ones(size(FH)); IBOUND(~isnan(FH))=-1;

% finite difference solution
[Phi,Q,Psi,Qx]=fdm2(xGr,yGr,kx,kx,IBOUND,FH,FQ);

%% show results
figure;

subplot(2,1,1);
hold on
plot(xm,FH( 1,:),'k');  % opgelegd
plot(xm,Phi(2:2:end,:),'b','linewidth',1);  % numeriek
plot(X, phi,'r')        % analytisch
plot(xPeilbuis,hPeilbuis,'ko','markerfacecolor','k');
grid on
title(sprintf('Verloop van de stijghoogte numeriek en analytisch, scenario %d',scenario));
legend('polderpeil','stijgh. numeriek','stijgh. analytisch','peibuis',4);
xlabel('x West-Oost [m] tov hart ARK');
ylabel('m +NAP');

%cnum=Dy(1)./kx(1,:)/2; cnum(cnum>1e5)=NaN;
%figure; plot(xm,cnum,'m'); legend('c');

%% check discharge in x-direction

%figure; hold on
subplot(2,1,2);
hold on;
plot(xGr(2:end-1),Qx(2,:),'b');  % numerical  soltuion
plot( X, q,'r'); grid on         % analytical solution
legend('Qx numerical','Qx analytical',4);

title(sprintf('Horizontal discharge through aquifer (eastward positive), scenario %d',scenario));
ylabel('q [m2/d]');
xlabel('x West-East [m] realtive of heart line Amsterdam Rhine Canal');

fprintf('Computed head in piezometere 31E0166 in %d is NAP %.2f m should be NAP %.2f m\n',...
          scenario,phi(find(X<=xPeilbuis,1,'last')),hPeilbuis);
      
pielbuis31E0163=phi(find(X<=1200,1,'last'));

