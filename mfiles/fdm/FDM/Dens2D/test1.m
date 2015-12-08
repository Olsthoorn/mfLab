% Test van PhiPsi.
% PhiPsi runt flatmeshctrd en psimeshctrd en plot the head and stream contours.
% Het kan zo alle interne bronnen aan, maar soms lijkt er toch een afwijking in het resultaat,
% namelijk dat de stijghoogtelijnen niet overal loodrecht staan op de stroomlijnen.
% dit probleem was op 000528 nog niet opgelost.
% TO000528

% De interne bronnen worden geinterpreteerd als vertical screens in de Psi berekening. Een onttrekking
% Q in het punt i,j is dan de lijnonttrekking direct onder dit punt. Dit heeft een vervelende consequentie
% voor een bovenrand met gegeven Q. Deze worden gezien als kleine screens. Hierdoor wijken de stroomlijnen
% in de bovenste rij cellen af. Dit probleem kan worden omzeild door de eerste cellenlaag zeer klein te kiezen.
% Niet zo elegant, maar het werkt. Het is vanwege de meerwaardigheid nu eenmaal nauwelijks mogelijk om
% eenduidig een stream function te berekenen voor een model met interne bronnen, zodanig dat die onder
% alle omstandigheden bevredigend is.
% TO000530

rdw=7.5;
x=[cosspace(0,rdw-0.1,20),rdw,dspace(rdw+0.1,10000,1.3,40)];
z=[-5.8,-6.2,-9.6,-9.7,-10.5,-10.6,-11.6,-11.7,-12,-12.3,-12.6,...
   cosspace(-12.8,-16,10),...
   cosspace(-16.1,-17.6,5),...
   cosspace(-17.7,-19.8,20),...
       -20,-25,-25.2,-200];
   
Nx=length(x);
Nz=length(z);

xm=(x(1:end-1)+x(2:end))/2;
zm=(z(1:end-1)+z(2:end))/2;

kh=20; kv=7;
ktop=0.01;
kEem=0.01;
kDW=1e-3;

kx=kh*ones(Nz-1,Nx-1);
I=find(zm<-4   ); kx(I,:)=ktop;  kz(I,:)=ktop;
I=find(zm<-9.6 ); kx(I,:)=10;    kx(I,:)=10;
I=find(zm<-10.6); kx(I,:)=ktop;  kz(I,:)=ktop;
I=find(zm<-11.6); kx(I,:)=kh;	  kz(I,:)=kv;
I=find(zm<-20  ); kx(I,:)=kEem;  kz(I,:)=kEem; 
I=find(zm<-25  ); kx(I,:)=kh;    kz(I,:)=kv;
I=find(zm>-16  ); J=find(xm>rdw & xm<rdw+0.1);  kx(I,J)=kDW; kz(I,J)=kDW;

FH=NaN*ones(Nz,Nx);
FH(1,:)=-6;
FH(find(zm<=-9.6 & zm>=-10.6),max(find(xm<rdw)))=-9.9;
FH(find(zm>=-17.6& zm<=-12.6),max(find(xm<rdw)))=-8.6;

FQ=zeros(size(FH));
%FQ(find(z>=-16 & z<=-11),find(abs(x-18)<1e-6))=10;

%FQ(74,30)=10;
NPhi=50; NPsi=50;
[Phi,Q,PsiStruct,PsiJump]=phipsi(x,z,kx,kz,FH,FQ,NPhi,NPsi);


set(gca,'xlim',[0,abs(min(z))]);