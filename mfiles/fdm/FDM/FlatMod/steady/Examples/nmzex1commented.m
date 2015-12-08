% PAO 2-3- nov 1999, NUMERIEK modelleren, T.N.Olsthoorn
% Voorbeeld: nmazex1 (nlagen mazure, example 1)
% Voorzien van printstatements en pauzeerinstructies

clc;				% clear screen
for i=1:30; fprintf('\n'); end
fprintf('Vergelijking van het nlagen mazure model met het 2d numerieke eindige differentiemodel\n');
fprintf('Drie-lagen systeem (dwz. 3 watervoerend pakketten met 3 scheidende lagen en boven vastgehouden\n');
kD=[1000;500; 750]
c =[ 400;800;1100]
Phi0=[-2;-1;-3]
pauzeer,

fprintf('Analytische oplossing n-lagen probleem:\n');
fprintf('Systeemmatrix: A=sysmat(kD,c); \n'); 
A=sysmat(kD,c)
pauzeer(5);

x=[0,logspace(0,4,50)]; dx=diff(x); Nx=length(x);

Phi=zeros(length(kD),length(x));
q  =zeros(length(kD),length(x));

fprintf('Berekenen Phi voor alle x-waarden en lagen');
for i=1:length(x)
   fprintf('.');
   Phi(:,i) =  expm(-x(i)*sqrtm(A))*Phi0;
   q  (:,i) = -diag(kD) * expm(-x(i)*sqrtm(A))* -sqrtm(A) *Phi0;   
end
fprintf('klaar\n');
pauzeer(5);

leg=[]; for i=1:size(Phi,1), leg=[leg,{['laag ',int2str(i)]}]; end
fprintf('Tonen berekende stijghoogten en fluxen in alle lagen\n');
subplot(2,1,1); plot(x,Phi); grid; title('head'); ylabel('[m]'); hold on; legend(leg);

subplot(2,1,2); plot(x,q);   grid; title('flow'); ylabel('[m2/d]'); xlabel('x [m]'); hold on; legend(leg);
pauzeer(5);

fprintf('Numerieke oplossing n-lagenprobleem:\n');
fprintf('z-waarden:\n');
dz=[5;20;20;35;10;50]; z=[0;-cumsum(dz)] %laagdikten en zwaarden van de knooppunten
Nz=length(z);
pauzeer(3);

fprintf('kwaarden berekend uit de kD en c waarden van de wvp''s en sdp''s:\n');
fprintf('k =[dz(1)/c(1);kD(1)/dz(2);dz(3)/c(2);kD(2)/dz(4);dz(5)/c(3);kD(3)/dz(6)]\n');
k =[dz(1)/c(1);kD(1)/dz(2);dz(3)/c(2);kD(2)/dz(4);dz(5)/c(3);kD(3)/dz(6)]
kx=k*ones(1,Nx-1);
kz=kx;

fprintf('Gegeven stijghoogten numerieke model, eerst matrix NaN''s dan fixed heads invullen linker en bovenrand:\n');
fprintf('FH=NaN*ones(Nz,Nx); FH([2,4,6],1)=Phi0; FH([3,5,7],1)=Phi0; FH(1,:)=0;\n');
FH=NaN*ones(Nz,Nx); FH([2,4,6],1)=Phi0; FH([3,5,7],1)=Phi0; FH(1,:)=0;
fprintf('Model runnen, mfile flatmeshctrd moet in het zoekpad van MATLAB staan of op de gebruikte directoty:\n');
fprintf('[Phi,Q]=flatmeshctrd(x,z,kx,kz,FH,[]);\n');
pauzeer(3);
tic
[Phi,Q]=flatmeshctrd(x,z,kx,kz,FH,[]);
fprintf('Het runnen kostte %f seconden\n',toc);
pauzeer(5);

fprintf('Numeriek berekende stijghoogten tekenen over het analytisch berekende verloop:\n');
fprintf('subplot(2,1,1); plot(x,Phi)\n');
subplot(2,1,1); plot(x,Phi,'+')
pauzeer(5),

fprintf('Uitstroming per laag berekenen voor x=0:\n');
Q0=[sum(Q(2:3,1));sum(Q(4:5,1));sum(Q(6:7,1))];  % flow per aquifer at x=0 numerical model
fprintf('Zet analytische en numerieke uitstroming naast elkaar:\n');
fprintf('   Analytic  Numeric\n');
disp([q(:,1),Q0]);									 % flow per aquifer at x=0 analytic model
pauzeer(3),


fprintf('\n\n\n\STROOMFUNCTIE:\n\n\n');
fprintf('Randvoorwaarden halen uit numeriek berekende Q:\n');
fprintf('Stroming van buiten het model in over de celwanden:\n');
Phiv=0.5*(Phi(:,1:end-1)+Phi(:,2:end));
fprintf('Instroming over de bovensrand:\n');
dqtop  =kz(1  ,:).*dx.*(Phiv(  1,:)-Phiv(    2,:))/dz(1);				% instroming over de bovenrand
fprintf('Idem over de onderrand:\n');
dqbot  =kz(end,:).*dx.*(Phiv(end,:)-Phiv(end-1,:))/dz(end);				% instroming over de onderrand

Phih=0.5*(Phi(1:end-1,:)+Phi(2:end,:));
fprintf('Instroming over de linker rand:\n');
dqleft =kx(  :,  1).*dz.*(Phih(:,  1)-Phih(:,    2))/dx(1);         % instroming over de linkerrand
fprintf('idem over de rechter rand:\n');
dqright=kx(  :,end).*dz.*(Phih(:,end)-Phih(:,end-1))/dx(end);       % instroming over de rechterrand


FPsi=NaN*ones(Nz,Nx);
FPsi(1      ,1      )=0;															% stel linker boven rand psi = 0.
FPsi(1      ,2:end  )=FPsi(  1,  1)+cumsum(dqtop);							% integratie Q langs bovenrand
FPsi(2:end  ,  end  )=FPsi(  1,end)+cumsum(dqright);							% idem langs rechterrand
FPsi(  end  ,1:end-1)=FPsi(end,end);											% idem langs dichte onderrand
FPsi(1:end-1,1      )=FPsi(end,  1)+flipud(cumsum(flipud(dqleft)));		% idem langs linkerrand
FPsi,
pauzeer(5);

fprintf('Berekening Psi met het gewone 2D model:\n');
pauzeer,
tic
[Psi,QPsi]=flatmeshctrd(x,z,1./kz, 1./kz,FPsi,[]);			% QPsi is alleen nuttig voor dichtheidsstroming, is hier dummy
fprintf('Berekenen stroomfunctie kostte %f seconden\n',toc);
pauzeer(5);

fprintf('Plotten van de stijghoogten en de stroomfunctie:\n');
fprintf('Accentueren scheidende lagen met patches:\n');
figure
patch([x(1),x(end),x(end),x(1)],[z(1),z(1),z(2),z(2)],'m'); hold on
patch([x(1),x(end),x(end),x(1)],[z(3),z(3),z(4),z(4)],'m');
patch([x(1),x(end),x(end),x(1)],[z(5),z(5),z(6),z(6)],'m');
pauzeer(3),

fprintf('Contouren van de stijghoogten:\n');
fprintf('contour(x,z,Phi,30);\n');
contour(x,z,Phi,30);
pauzeer(5),

fprintf('idem van de stroomlijnen:\n');
fprintf('c=contour(x,z,Psi,30);\n');
c=contour(x,z,Psi,30);
xlabel('x [m]'); ylabel('z [m]'); title('Doorsnede met 30 stijghoogte- en 30 stroomlijnen');
pauzeer(5),

fprintf('\n\n\nVERBLIJFTIJDEN:\n\n\n');
fprintf('Functie getcontours creert een contour-struct met fields Psivalue,N,x,z,phi,t:\n');
fprintf('C=getcontours(x,z,Phi,kx,kz,por,c);\n');

por=0.35*ones(size(kx));
tic
C=getcontours(x,z,Phi,kx,kz,por,c);
fprintf('Het berekenen van de verblijftijden kostte %f seconden.\n',toc);
pauzeer(5);

fprintf('Toon fields van contourstruct\n');
C
pauzeer(3);

fprintf('Zet alle tijden in een vector:\n');
T=cat(1,C.T);
fprintf('Maak grafiek van de cumulatieve verblijftijd van alle stroomlijnen:\n');
fprintf('figure\n');
fprintf('plot(sort(T),[1:length(T)]/length(T)*100);\n');
figure
plot(sort(T),[1:length(T)]/length(T)*100);
xlabel('verblijftijd [d]'); ylabel('% volumestroom'); title('cumulatieve verblijftijd');
