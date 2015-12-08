% PAO 2-3- nov 1999, NUMERIEK modelleren, T.N.Olsthoorn
% Voorbeeld: nmazex1 (nlagen mazure, example 1)

clc;				% clear screen
for i=1:30; fprintf('\n'); end
%Vergelijking van het nlagen mazure model met het 2d numerieke eindige differentiemodel\n');
%Drie-lagen systeem (dwz. 3 watervoerend pakketten met 3 scheidende lagen en boven vastgehouden\n');
kD=[1000;500; 750]
c =[ 400;800;1100]
Phi0=[-2;-1;-3]

%Systeemmatrix: A=sysmat(kD,c); \n'); 
A=sysmat(kD,c)

x=[0,logspace(0,4,50)]; dx=diff(x); Nx=length(x);

Phi=zeros(length(kD),length(x));
q  =zeros(length(kD),length(x));

%Berekenen Phi voor alle x-waarden en lagen');
for i=1:length(x)
   Phi(:,i) =  expm(-x(i)*sqrtm(A))*Phi0;
   q  (:,i) = -diag(kD) * expm(-x(i)*sqrtm(A))* -sqrtm(A) *Phi0;   
end

leg=[]; for i=1:size(Phi,1), leg=[leg,{['laag ',int2str(i)]}]; end
subplot(2,1,1); plot(x,Phi); grid; title('head'); ylabel('[m]'); hold on; legend(leg);
subplot(2,1,2); plot(x,q);   grid; title('flow'); ylabel('[m2/d]'); xlabel('x [m]'); hold on; legend(leg);

dz=[5;20;20;35;10;50]; z=[0;-cumsum(dz)]; %laagdikten en zwaarden van de knooppunten
Nz=length(z);

k =[dz(1)/c(1);kD(1)/dz(2);dz(3)/c(2);kD(2)/dz(4);dz(5)/c(3);kD(3)/dz(6)]
kx=k*ones(1,Nx-1);
kz=kx;

FH=NaN*ones(Nz,Nx); FH([2,4,6],1)=Phi0; FH([3,5,7],1)=Phi0; FH(1,:)=0;	% fixed heads model
[Phi,Q]=flatmeshctrd(x,z,kx,kz,FH,[]);			% run model fixed Q's just [ ]

subplot(2,1,1); plot(x,Phi,'+')

Q0=[sum(Q(2:3,1));sum(Q(4:5,1));sum(Q(6:7,1))];  % flow per aquifer at x=0 numerical model
fprintf('   Analytic  Numeric\n');
disp([q(:,1),Q0]);									 % flow per aquifer at x=0 analytic model


%Stroomfunctie:\n');
%Stroming van buiten het model in over de celwanden:\n');
Phiv=0.5*(Phi(:,1:end-1)+Phi(:,2:end));
dqtop  =kz(1  ,:).*dx.*(Phiv(  1,:)-Phiv(    2,:))/dz(1);				% instroming over de bovenrand
dqbot  =kz(end,:).*dx.*(Phiv(end,:)-Phiv(end-1,:))/dz(end);				% instroming over de onderrand

Phih=0.5*(Phi(1:end-1,:)+Phi(2:end,:));
dqleft =kx(  :,  1).*dz.*(Phih(:,  1)-Phih(:,    2))/dx(1);         % instroming over de linkerrand
dqright=kx(  :,end).*dz.*(Phih(:,end)-Phih(:,end-1))/dx(end);       % instroming over de rechterrand

FPsi=NaN*ones(Nz,Nx);
FPsi(1      ,1      )=0;															% stel linker boven rand psi = 0.
FPsi(1      ,2:end  )=FPsi(  1,  1)+cumsum(dqtop);							% integratie Q langs bovenrand
FPsi(2:end  ,  end  )=FPsi(  1,end)+cumsum(dqright);							% idem langs rechterrand
FPsi(  end  ,1:end-1)=FPsi(end,end);											% idem langs dichte onderrand
FPsi(1:end-1,1      )=FPsi(end,  1)+flipud(cumsum(flipud(dqleft)));		% idem langs linkerrand

[Psi,QPsi]=flatmeshctrd(x,z,1./kz, 1./kz,FPsi,[]);			% QPsi is alleen nuttig voor dichtheidsstroming, is hier dummy

%Plotten van de stijghoogten en de stroomfunctie:\n');
figure
patch([x(1),x(end),x(end),x(1)],[z(1),z(1),z(2),z(2)],'m'); hold on
patch([x(1),x(end),x(end),x(1)],[z(3),z(3),z(4),z(4)],'m');
patch([x(1),x(end),x(end),x(1)],[z(5),z(5),z(6),z(6)],'m');
  contour(x,z,Phi,30);
c=contour(x,z,Psi,30);
xlabel('x [m]'); ylabel('z [m]'); title('Doorsnede met 30 stijghoogte- en 30 stroomlijnen');

%Verblijftijden =====================================================================\n');
por=0.35*ones(size(kx));
C=getcontours(x,z,Phi,kx,kz,por,c);
T=cat(1,C.T);
figure
plot(sort(T),[1:length(T)]/length(T)*100);
xlabel('verblijftijd [d]'); ylabel('% volumestroom'); title('cumulatieve verblijftijd');



