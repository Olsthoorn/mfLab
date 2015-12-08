% Papendorp bemaling sleuf 40 m lang op 30 m van het ARK kanaal
% veel laagjes in de verticaal
%
% TO 001005
close all
clear

cA=40;		% c bodem ARk-kanaal
cS=0.01		% c bodem sleuf
cP=500;		% c polder, deklaag
kh=45;		% WVP
kv=45000;		% WVP
k=sqrt(kh*kv); % mix
z=[0;-3.5;cosspace(-3.75,-53.5,20,'R')']; dz=abs(diff(z)); zm=0.5*(z(1:end-1)+z(2:end));

% model 1
x=	[fliplr(-logspace(log10(140),log10(4000),10)),cosspace(-200,-140,10,'L'),...
      cosspace(-140,0,20),...
      cosspace(0,30,15),...
      cosspace(30,33,7),...
      cosspace(33,100,15,'R'),logspace(log10(33),log10(4000),10)];
x=unique(round(100*x))/100;
xm=0.5*(x(1:end-1)+x(2:end));

Nx=length(x)-1; Nz=length(z)-1;			% Nx,Nz is number of cells!

% find cells of ARK en sleuf
IAm=find(xm>-140 & xm<0);
ISm=find(xm>30 & xm<33);

% find cells of layers
JSDL=find(zm<z(1) & zm>z(2));		% aquitard
JWVP=find(zm<z(2));					% aquifer

KH=kh*ones(Nz,Nx);
KH(JSDL,:)=1e-3;						% k van polder toplaag
KH(JWVP,:)=kh;							% k wvp

KV=KH;
KV(JSDL,:)=0.5*dz(1)/cP;					% k van polder toplaag
KV(1,IAm) =0.5*dz(1)/cA;					% k onder ARK
KV(1,ISm) =0.5*dz(1)/cS;					% k sleuf

% fixed heads
FH=NaN*zeros(Nz,Nx);
FH(1,:)=0;
FH(1,ISm)=-2.1;
FQ=zeros(size(FH));

[Phi,Q]=flatblockctrd(x,z,KH,KV,FH,FQ);
xlim=[-200,100];

subplot(3,1,1); plot(xm,Phi,'-xb');							set(gca,'xlim',xlim); grid on; title('volkomen, phi');
subplot(3,1,2); plot(x,[0,cumsum(Q(1,:))],'-xb');		set(gca,'xlim',xlim); grid on; title('volkomen, q');
subplot(3,1,3); plot(xm,Q(1,:)./diff(x),  '-xr');		set(gca,'xlim',xlim); grid on; title('volkomen, v');
qA=sum(Q(1,IAm),2);
qS=sum(Q(1,ISm),2);

figure; contour(xm,zm,Phi,30); set(gca,'xlim',[-200,100]);

disp([qA,qS]);














xx=[-140,0,30,33];					% sectie coordinaten
h=[0,0,0,-2.1,0];						% polderpeil, kanaalpeil, polderpeil, bouwput(bouwsleufpeil), en polderpeil
c=[500];									% weerstand deklaag  (we hebben maar een scheidende laag, de deklaag)
T=[41*55];								% weerstand WVP (we hebben maar een wvp)
z=[-0.5,-3.5,-55]';					% vlakhoogten voor het tekenen, let op de " ' " die er een vertical kolom van maakt
											% De vlakhoogten zijn: top deklaag, onderzijde deklaag onderzijde WVP
c=[c,c,c,c,c];							% elke sectie zelfde c -waarde
T=[T,T,T,T,T];							% elke sectie zelfde kD-waarde
c(1,2)=40; c(1,4)=0.0001;			% c-waarde (1,2) wordt kanaalc-waarde, c(1,4) wordt dit van de sleuf (ongeveer nul)

Q=zeros(size(c(:,1:end-1)));		% geen onttrekkingen op de grensvlakken van de secties

X=[-200:1:100];						% tekenen tusse x=-300 en x=300 m

[Section,phi,q,s,X,xx,kD,c,h,Q]=nsecn(xx,T,c,h,0*h,Q,X);			% analytische oplossing voor alle X waarden en lagen

shownsecn(phi,q,s,z(:,1),h,X,xx)							% laat maar zien
set(gca,'xlim',[X(1),X(end)],'ylim',[-55,0]);		% pas assen aan aan deze situatie

%Met weerstand sleuf
xx=[  -140,  0,  30,  33];			% sectie coordinaten
h=[0,0,0,-2.1,0];						% polderpeil, kanaalpeil, polderpeil, bouwput(bouwsleufpeil), en polderpeil
c=[500];									% weerstand deklaag  (we hebben maar een scheidende laag, de deklaag)
T=[41*55];								% weerstand WVP (we hebben maar een wvp)
z=[-0.5,-3.5,-55]';					% vlakhoogten voor het tekenen, let op de " ' " die er een vertical kolom van maakt
											% De vlakhoogten zijn: top deklaag, onderzijde deklaag onderzijde WVP
c=[c,c,c,c,c];							% elke sectie zelfde c -waarde
T=[T,T,T,T,T];							% elke sectie zelfde kD-waarde
c(1,2)=40; c(1,4)=0.062;			% c-waarde (1,2) wordt kanaalc-waarde, c(1,4) wordt dit van de sleuf (ongeveer nul)

Q=zeros(size(c(:,1:end-1)));		% geen onttrekkingen op de grensvlakken van de secties

X=[-200:1:100];						% tekenen tusse x=-300 en x=300 m
[Section,phi1,q1,s1,X,xx,kD,c,h,Q]=nsecn(xx,T,c,h,0*h,Q,X);	% analytische oplossing voor alle X waarden en lagen

shownsecn(phi1,q1,s1,z(:,1),h,X,xx)							% laat maar zien
set(gca,'xlim',[X(1),X(end)],'ylim',[-55,0]);		% pas assen aan aan deze situatie
