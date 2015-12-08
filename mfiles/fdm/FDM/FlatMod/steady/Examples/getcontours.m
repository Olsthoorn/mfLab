function C=getcontours(X,Z,Phi,kx,kz,Por,c)
% puts contours c in a struct capital C containing value, number of coordinates, x and y values

i=1; i2=0;

% uiteenrafelen van de contouren, afzonderlijke contourlijnen worden in de struct C geplaatst
while 1
i1=i2+1;
C(i).value=c(1,i1);
C(i).N    =c(2,i1);
i1=i2+2;
i2=i2+1+C(i).N;
C(i).x=c(1,i1:i2);
C(i).z=c(2,i1:i2);
i=i+1;
if i2==length(c), break; end
end

% toevoegen middelpunten van de lijnstukjes van de contouren
% van de doorlatendheid en de porositeit op deze middenpunten
% van de gradienten van de lijnstukjes en de verblijftijd
Nk=length(kx(:,1));
for i=1:length(C)
   x=C(i).x;    z=C(i).z;   ds=sqrt(diff(x).^2+diff(z).^2);
   C(i).phi=interp2(X,Z,Phi,x,z);
   xM=0.5*(x(1:end-1)+x(2:end));
   zM=0.5*(z(1:end-1)+z(2:end));
   [I,J]=getcell(X,Z,xM,zM);   IJ=(J-1)*Nk+I;
   kh=kx(IJ);      kv=kz(IJ);   por=Por(IJ);
   dt=-ds./(diff(C(i).phi)./por.*kh);				% should be some combination of kh and kv, too difficult for now
   C(i).t=[0,cumsum(dt)];
   if min(dt)<0, C(i)=flipline(C(i)); end
   C(i).T=max(C(i).t);
end

function [I,J]=getcell(X,Z,xM,zM);
% selecteert de cell in netwerk bepaald door X,Z, waarin punt xM,zM ligt
I=zeros(size(xM));
J=zeros(size(zM));
for i=1:length(xM)
     I(i)=sum(X<=xM(i));
     J(i)=sum(Z>=zM(i));
end
  
function C=flipline(C)
% flips de streamline, zodat hij altijd met de tijd meeloopt.
C.x  = fliplr(C.x);
C.z  = fliplr(C.z);
C.phi= fliplr(C.phi);
C.t  = fliplr(C.t)-min(C.t);