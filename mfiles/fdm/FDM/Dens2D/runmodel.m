% runmodel
% Flat model to be used for flowanalyse (slem)
% TO 990520

Figdir='z:\projects\hy171 zoet zout modellering\psiphimodel\newfigs\';

%profielLeiduin;								% data for profile CC, Leiduin
profielDeZilkNew;								% data for profile EE, De Zilk

% The numerical model setup
xMin    =-5000;  xMax=15000;				% x-range limit
xZone(1)=xMin; xZone(end)=xMax;			% cut off the x-zones
DZ      =5;										% vertical cell size of numerical model
Y       =[yZone(1):-DZ:yZone(end)]';	% create numerical model layers
X       =xMin:25:xMax;						% X boundaries of numerical cells
dx      =diff(X);								% Numerical model cell widths
Nx      =length(X);							% total number of x-points in numerical model
Ny      =length(Y);							% tolal number of y-points in numerical model

if Y(end)>Y(1),								% make sure top row is always top of model
   Y=flipud(Y);
end

xc      =(X(1:end-1)+X(2:end))/2;		% cell center
yc      =(Y(1:end-1)+Y(2:end))/2;		% cell center

% in what hydrological zone lies each numerical cell center?
IxSect=zeros(size(xc)); for i=1:length(xc),   IxSect(i)=find(xc(i)>xZone(1:end-1) & xc(i)<xZone(2:end)); end
IySect=zeros(size(yc)); for i=1:length(yc),   IySect(i)=find(yc(i)<yZone(1:end-1) & yc(i)>yZone(2:end)); end

%we use precipitation everywhere and fixed head to set the heads at the top row of the model
%where we take the c1 as the entry resistance

QtopCell=dx.*prev(IxSect);													% Precipitation per dx
FQ=zeros(Ny,Nx);																% complete matrix of fixed Q's
FQ(1,:)=0.5*([0,QtopCell]+[QtopCell,0]);								% Put preciptation in top nodes
for j=1:size(QZ,1)															% for all aquifers
   for i=1:size(QZ,2)														% for all sections
      if QZ(j,i)
	      ix=find(X==x(i));													% find X node of section boundary
         iy=find(Y==yZone(2*j));											% find Y node of aquifer
         switch j																% switch aquifer
         case 1,																% top aquifer
            FQ(iy-1,ix)=FQ(iy-1,ix)+QZ(j,i);							% phreatic extraction into top node
            zWF=[0;-5];														% extraction depth
         case 2,																% second aquifer, well screen length
            FQ(iy+1,ix)=FQ(iy+1,ix)+QZ(j,i)*(1/2);					% z= -25 m
            FQ(iy+3,ix)=FQ(iy+3,ix)+QZ(j,i)*(1/2);					% z= -35 m
            zWD=[-20;-40];													% extraction depth
         otherwise
            FQ(iy,ix)=FQ(iy,ix)+QZ(j,i);								% all other (=third) aquifers
         end
      end
   end
end

% fixed heads at the top of the system
FPhi=zeros(1,Nx)*NaN;														% start with the top nodes
H =h(IxSect);																	% H of the sections, not the nodes
I=1+find(~isnan(H(1:end-1))&~isnan(H(2:end))); FPhi(I)=0.5*(H(I-1)+H(I));	% left and right fixed heads, no nans
I=1+find( isnan(H(1:end-1))&~isnan(H(2:end))); FPhi(I)=            H(I);	% left NaN, right noNaN, use right
I=1+find(~isnan(H(1:end-1))& isnan(H(2:end))); FPhi(I)=     H(I-1);			% right NaN, left noNaN, use left
FPhi(1)=H(1); FPhi(end)=H(end);											% first and last node
FPhi=[FPhi;NaN*zeros(Ny-1,Nx)];											% expand the matrix to all nodes

% entries at the top of the system
CE=zeros(1,Nx)*NaN;															% start with the top nodes
ce=dx./cE(IxSect);															% entry conductance per section, not node
I=1+find(~isnan(ce(1:end-1))&~isnan(ce(2:end))); CE(I)=0.5*(ce(I-1)+ce(I)); % left and right have entry
I=1+find( isnan(ce(1:end-1))&~isnan(ce(2:end))); CE(I)=0.5*         ce(I) ; % left NaN, right noNan
I=1+find(~isnan(ce(1:end-1))& isnan(ce(2:end))); CE(I)=0.5* ce(I-1);			 % right NaN, left noNaN
CE(1)=ce(1); CE(end)=ce(end);												% first and last
CE=[CE;NaN*zeros(Ny-1,Nx)];												% expand the matrix to all nodes


%Other specified fixed heads, elsewhere in the model
for j=1:size(FH,1)															% for all aquifers
   for i=1:size(FH,2)														% for all sections
      if ~isnan(FH(j,i))													% if not a NaN
	      ix=find(X==x(i));													% find X node of section boundary
         iy=find(Y==yZone(2*j));											% find Y node of aquifer
         switch j																% switch aquifer
         case 1,																% top aquifer
            FPhi(iy-1,ix)=FH(j,i);										% phreatic extraction into top node
            CE(  iy-1,ix)=NaN;
            zWF=[0;-5];														% extraction depth
         case 2,																% second aquifer, well screen length
            FPhi(iy+1,ix)=fh(iy+1,ix);									% z= -25 m
            FPhi(iy+2,ix)=fh(iy+2,ix);									% z= -30 m
            FPhi(iy+3,ix)=fh(iy+3,ix);									% z= -35 m
         otherwise
            FPhi(iy,ix)=fh(iy,ix);											% all other (=third) aquifers
         end
      end
   end
end

if bufzone				% bufferzone
i=find(X==4225); j=[1:4]'; FPhi(j,i)=phi1994(j,i); CE(j,i)=NaN;
end

SMALL=1e-3;
kx      =kLayers(IySect,IxSect);			% kx
ky      =kLayers(IySect,IxSect);  		% ky=kx
kx(1,find(xc>x(end)))=SMALL;				% eliminate first model cells in Haarlemmermeer

tic
DX=125; DY=5;												% horizontal and vertical interface thickness (erfc argument scaler)
fprintf('computing density term');
 [dddx,dddy]=densterm(X,Y,isodens,ky,DX,DY);		% densterm  uses erfc function for interface
%[dddx,dddy]=densterm1(X,Y,isodens,ky);			% densterm1 uses sharp interfaces
toc

tic
fprintf('running the phimodel...\n');
[Phi,Q]=flatmeshctrd1(X,Y,kx,ky,FPhi,FQ+dddx,CE);			   % Ready, including density effect
QPhi=Q-dddx;
fprintf('running the psimodel...\n');
[Psi,QPsi,FPsi,Qcum]=psimeshctrd(X,Y,kx,ky,QPhi,dddy);

% find internal sources:
iSS=[]; mustsplit=any(QZ~=0|~isnan(FH));
for i=1:length(x)
   if mustsplit(i) | (bufzone & x(i)==4225),						% bufzone is switch for bufzone, 4225m=dune rim
      iSS=[iSS,find(X==x(i))];
   end 
end
iSS=[iSS,length(X)];

clear PsiStruct
i1=1;
for i=1:length(iSS);
   i2=iSS(i);
   PsiStruct(i).Psi=Psi(:,i1:i2);
   PsiStruct(i).Psi(:,1)=PsiStruct(i).Psi(:,1)-Qcum(:,i1);
   PsiStruct(i).x=X(i1:i2);
   i1=i2;
end
toc

tic
fprintf('plotting results...\n');
FS=12;

figure; axes; hold on                         % Hydrosomes of brackish and saline water first
yLim=[yZone(end),35];

[xiso1,yiso1]=clip(X,Y,[isodens(1).x;isodens(1).x(end);isodens(1).x(1)]',[isodens(1).y;min(Y);min(Y)]');
	fill(xiso1,yiso1,[0 1 0]);		% Brackish hydrosome
hold on
[xiso2,yiso2]=clip(X,Y,[isodens(2).x;isodens(2).x(end);isodens(2).x(2)]',[isodens(2).y;min(Y);min(Y)]');
	fill(xiso2,yiso2,[0 1 1]);		% Salt hydrosome
   
%fill([X(1),X(end),X(end),X(1)],[yLim(2),yLim(2),0,0],'w');			% cover top of hydrosomes   
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'hydrsm',int2str(T),'.png''']);
close(gcf);

%plot layout, 
if ir==1
figure; axes; hold on
plot(Duin(:,end),Duin(:,end-1),'k');							% dune profile

%coloring c-zones
cLayers=[1,3,5];
for iL=cLayers; dd=0;
   for iz=1:length(h)
      if ~(iL==1 &  (isnan(h(iz)) | xZone(iz)<0))
         if iL==1 & xZone(iz+1)>xZone(end-1), dd=-5; else dd=0; end
         fill([xZone(iz),xZone(iz+1),xZone(iz+1),xZone(iz)],...
              [yZone(iL),yZone(iL),yZone(iL+1),yZone(iL+1)]+dd,...
               zoneColor(kLayers(iL,iz),kLayers(cLayers,:)));
       end
    end
 end
 for iz=1:length(x);
    if QZ(1,iz)|~isnan(FH(1,iz)),
       plot([x(iz),x(iz)],[zWF(1),zWF(2)],'b');				% plot phreatic extraction
    end
    if QZ(2,iz)|~isnan(FH(2,iz)),
       plot([x(iz),x(iz)],[zWD(1),zWD(2)],'b');				% plot extraction wells
    end
 end
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'hydrsm',int2str(T),'.png''']);
close(gcf);

% afdekplaatje op haarlemmermeer
figure; axes; hold on
iL=1; iz=length(h);
fill([xZone(iz),xZone(iz+1),xZone(iz+1),xZone(iz)],...
     [yZone(iL),yZone(iL),yZone(iL+1),yZone(iL+1)],'m');
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'coverHmeer',int2str(T),'.png''']);
close(gcf);
end



% contouring phi
figure; axes; hold on
if 1
	fprintf('contouring Phif\n');
   dphi=0.25;																% phi increment
   PhiRange=[floor(min(Phi(:))):dphi:ceil(max(Phi(:)))];		% Range to plot contours
   contour(X,Y,Phi,PhiRange,'r');									% Isohypses
end
if 0
   rhos=1.022; rhof=1.00;
   fprintf('contouring PhiS\n');
   Phis=rhof/rhos*(Phi-Y*ones(size(X)))+Y*ones(size(X));
   dphi=0.25;																% phi increment
   PhiSRange=[floor(min(Phis(:))):dphi:ceil(max(Phis(:)))];	% Range to plot contours
   contour(X,Y,Phis,PhiSRange,'r');									% Isohypses
end
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'philines',int2str(T),'.png''']);
close(gcf);

fprintf('Contouring psi\n');
figure; axes; hold on
for i=1:length(PsiStruct)
   dpsi=0.10;																% Psi increment (m2/d)
   PsiRange=[floor(min(PsiStruct(i).Psi(:))):dpsi:ceil(max(PsiStruct(i).Psi(:)))];
   contour(PsiStruct(i).x,Y,PsiStruct(i).Psi,PsiRange,'k');	% IsoPsi lines
   hold on
end
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'psilines',int2str(T),'.png''']);
close(gcf);

fprintf('heads in cross section touch\n');
figure; axes; hold on;
   plot(X,Phi(1,:),'k');												% phreatic
   plot(X,Phi(3,:),'k');												% first aquifer
   plot(X,Phi(5,:),'r');												% second aquifer
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
eval(['print -dpng ''',Figdir,'phi1_3_5',int2str(T),'.png''']);
close(gcf);

figure; axes; hold on
set(gca,'yLim',yLim,'xlim',[X(1),X(end)],'fontsize',FS);
title(sprintf('EW-section %s AWD, dPhi=%4.2fm, dPsi=%4.2fm2/d (yr=%d)',ProfileName,dphi,dpsi,T),'fontsize',FS); 
xlabel('Distance from shore [m]','fontsize',FS);
ylabel('Elevation [m +msl]','fontsize',FS);
eval(['print -dpng ''',Figdir,'title',int2str(T),'.png''']);
close(gcf);

toc   

if 0
	figure
	J=find(X<=9000 & X>=-1000);
	%plot(Duin(:,end),Duin(:,end-1),'k'); hold on
	plot(X(J),Phi(find(Y== -5),J),'b',...
   X(J),Phi(find(Y==-20),J),'r',...
   X(J),Phi(find(Y==-75),J),'g');
	grid
	legend('z=- 5m','z=-20m','z=-75m');
	xlabel('x from coast [m]'); ylabel('elevation [m +NAP]'); title(['heads (',int2str(T),')']);
   eval(['save phi',int2str(T),' Phi']);
   vls=Phi(1,find(X==x(2)))
	okf=Phi(1,find(X==x(5)))
	okd=Phi(find(Y==-30),find(X==x(5)))
	drd=Phi(find(Y==-30),find(X==x(6)))
end
toc