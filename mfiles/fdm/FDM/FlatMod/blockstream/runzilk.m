% test model blkstrm, block centred finite difference model with stream function
% TO 000530, 001026, stroomfunctie ingebouwd.

profielDeZilk;

X=unique([-5000:25:8000,x]);	DX=diff(X,1,2);
xZone=[X(1),x,X(end)];
yZone =[0;-cumsum(D)];						% vertical analytic layer boundaries, elevations
Y=[yZone]'; Y=Y(:);					% leg laag middens in midden en op dranden van wvp's

Nx=length(X);
Ny=length(Y);

Xm=0.5*(X(2:end)+X(1:end-1));
Ym=0.5*(Y(2:end)+Y(1:end-1));

izone=zeros(size(Xm));
for i=1:length(xZone)-1,    izone(find(Xm>xZone(i) & Xm<xZone(i+1)))=i;		end

k=kLayers(:,izone); 

FH=zeros(Ny-1,Nx-1)*NaN; FH(1,:)=h(izone);
FQ=zeros(Ny-1,Nx-1);		FQ(1,:)=prev(izone).*DX(1,:);

[phi,Q,Psi]=blckstrm(X,Y,k,k,FH,FQ);
Phi=zeros(length(yZone),length(Xm));
Phi(1:Ny-1,:)=phi;
Phi(3:2:Ny,:)=phi(2:2:end,:)

close all
figure; axes; hold on
for i=1:length(xZone)-1
   patch([xZone(i),xZone(i+1),xZone(i+1),xZone(i)],[yZone(4),yZone(4),yZone(3),yZone(3)],[0.8,0.8,0.8]);
   patch([xZone(i),xZone(i+1),xZone(i+1),xZone(i)],[yZone(6),yZone(6),yZone(5),yZone(5)],[0.8,0.8,0.8]);
end
contour(Xm,Y,Phi,50,'b')
contour(X(2:end-1),Y,Psi(:,2:end-1),200,'r');
plot(EE(:,4),EE(:,3));
