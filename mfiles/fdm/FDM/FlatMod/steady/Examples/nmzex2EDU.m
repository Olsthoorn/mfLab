% nmazex2 (nlagen mazure, example 1)
% PAO 2-3 november 1999, NUMERIEK MODELLEREN, mazure voorbeeld 2
%
% Situatie van Wassen (1990), doorsnede vanaf de Utrechtse Heuvelrug door de Bethunepolder
% tot in de Loosdrechtse plassen
% Theo Olsthoorn, 991014

%Het netwerk

x=[-2500	-1800	-1400	-1200	-1000	-900	-700	-500	-300	-100,...
   100	300	500	700	900	1000	1200	1400	1600	1800,...
      2200	2600	3000	3400	3800	4200	4600	5000	5600	6000,...
      6400	6800	7200	7600	8000	8400	8800	9200	9600	10000	10400	11000];

z=[0;-5;-10;-15;-25;-35;-40;-45;-50;-55;-60;-65;-70;-75;-80;-90;-100;-110;-130;-150];

dx=diff(x);    Nx=length(x);
dz=-diff(z);   Nz=length(z);

%Bovenranvoorwaarde, afgelezen uit fig. 2, p84 van Wassen (1990)
h=NaN*zeros(size(x));
I=find(x<-1000);            h(I)=-1.10;
I=find(x>-1000 & x<1000);   h(I)=-3.85;
I=find(x> 1000 & x<3250);   h(I)=-1.20;
I=find(x> 3250 & x<4500);   h(I)=-1.00;
I=find(x> 4500 & x<5500);   h(I)=-0.80;
I=find(x> 5500 & x<6500);   h(I)=-0.40;
I=find(x> 6500 & x<7250);   h(I)=-0.00;
I=find(x> 7250 & x<8750);   h(I)= 0.40;
I=find(x> 8750 & x<9750);   h(I)= 0.80;
I=find(x> 9750 & x< 10500); h(I)= 1.20;
I=find(x>10500 & x<=11000); h(I)= 1.60;
   
load 'idkBethuneEDU.txt'										% lees matrix met k-indices.
K     =[30;  4/10;  30/400;  30/400;  1/50];			% k-waarden bij id nummers in id matrix
Map   =[0.9; 0.8; 0.6; 0.6; 0.85]*[1,1,1];         % colormap voor inkleuring cellen naar k-waarden

k=zeros(Nz-1,Nx-1);
k=reshape(K(idkBethuneEDU(:)),Nz-1,Nx-1);					% kwaarden aanmaken met indices en K-lijstje
      
% gegeven stijghoogten numerieke model, eerst matrix NaN's dan fixed heads invullen linker en bovenrand
FH=NaN*ones(Nz,Nx); FH(1,:)=h;

% model runnen, mfile flatmeshctrd staat bij mij in subdirectory "flatmodel"
% tic, toc meet de rekentijd
tic
[Phi,Q]=flatmeshctrd(x,z,k,k,FH,[]);
fprintf('Runnen van het Phi-model (%d x %d knopen) kostte %f seconden\n',Nz,Nx,toc);

% stream function
% randvoorwaarden halen uit numeriek berekende Q
% stroming van buiten het model in over de celwanden:
Phiv=0.5*(Phi(:,1:end-1)+Phi(:,2:end));
dqtop  =k(1  ,:).*dx.*(Phiv(  1,:)-Phiv(    2,:))/dz(1);				% instroming over de bovenrand


FPsi=NaN*ones(Nz,Nx);
FPsi(1      ,1      )=0;																% stel linker boven rand psi = 0.
FPsi(1      ,2:end  )=FPsi(  1,  1)+cumsum(dqtop);							% integratie Q langs bovenrand
FPsi(2:end  ,  end  )=FPsi(  1,end);											% rechter rand gesloten
FPsi(  end  ,1:end-1)=FPsi(end,end);											% linker rand gesloten
FPsi(1:end-1,1      )=FPsi(end,  1);											% onderrand gesloten

% berekening Psi met het gewone 2D model
tic
[Psi,QPsi]=flatmeshctrd(x,z,1./k, 1./k,FPsi,[]);			% QPsi is alleen nuttig voor dichtheidsstroming, is hier dummy
fprintf('Runnen van het Psi-model (%d x %d knopen) kostte %f seconden\n',Nz,Nx,toc);

% plotten van de stijghoogten en de stroomfunctie
figure
colormap(Map);												% gebruik boven opgegeven kleurenschema
contour(x,z,Phi,[-4:0.1:4]);							% eerst contouren anders zet image plaatje op zijn kop (waarom??)
hold on
image(x,z,idkBethuneEDU);									% inkleuren modelcellen met k-waarden
contour(x,z,Phi,[-4:0.1:4],'r');						% rode stijghoogtecontouren tussen NAP -4 en +4 met stappen van 0.1
c=contour(x,z,Psi,60,'b');								% 60 blauwe stroomlijnen 
plot(x,h)													% stijghoogte verloop bovenrand als opgegeven
xlabel('x [m]'); ylabel('z [m]'); title('Doorsnede Heuvelrug, Bethune Loosdrechtse Plassen');

% verblijftijden, functie getcontours creert een contour-struct met fields Psivalue,N,x,z,phi,t.
por=0.35*ones(size(k));
C=getcontours(x,z,Phi,k,k,por,c);					% zie file nmzex1
T=cat(1,C.T)/365.25;										% haal verblijtijd uit C en zet om naar jaren

figure														% plot cumulatieve verblijftijdsverdeling
plot(sort(T),[1:length(T)]/length(T)*100);
xlabel('verblijftijd [jaren]');
ylabel('% totale volumestroom');
title('cumulatieve verblijftijd, [over alle stroomlijnen]');



