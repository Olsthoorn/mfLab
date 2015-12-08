% VES simulation using radmod
% TO 060127 for CT4420
close all
D=1;
Lmin=2*D; Lmax=250;  llmin=log10(Lmin); llmax=log10(Lmax);
rMin=0.1; rMax=1000; lrmin=log10(rMin); lrmax=log10(rMax);
zMin=0; zMax=1000;

L=logspace(llmin,llmax,31);
%r=unique(sort([logspace(lrmin,lrmax,41), L-D, L+D]));
r=logspace(-1,3,51);

% layer properties
LayTop=  [ 0 5 10 20 40 60 100 200];
resis =  [ 1 10 10 0.2 0.2 0.2 0.2 0.2]; 
Z=4*LayTop(end);

% build input for mdoel
z=unique(sort([linspace(0,Z,51), LayTop]))';
%z=unique(sort([0, logspace(-1,log10(zMax),61), LayTop]))';

zm=0.5*(z(1:end-1)+z(2:end));

k=zeros(length(z)-1,length(r)-1);
for i = 1:length(LayTop)
    k(find(zm>LayTop(i)),:)=1./resis(i);
end

FQ=zeros(length(z),length(r)); FQ(1,1)=1;
FH=NaN*FQ; FH(end,:)=0; FH(:,end)=0;

[Phi,Q]=radmod(r,z,k,k,FH,FQ);


%xVes=[-fliplr(r),r(1:end)]; zVes=z; PhiVes=[fliplr(Phi),Phi];

%PV=interp2(xVes+L(end),zVes,PhiVes,xVes,zVes)-interp2(xVes-L(end),zVes,PhiVes,xVes,zVes);
%figure
%hold on
%for i=1:length(zVes)
%    line([xVes(1),xVes(end)],-[zVes(i),zVes(i)]);
%end
%for i=1:length(xVes)
%    line([xVes(i),xVes(i)],-[zVes(1),zVes(end)]);
%end
%contour(xVes,-zVes,PV,20);

DV=[]; rho=[];
for i=1:length(L)
    r1=L(i)-D; r2=L(i)+D;
    DV(i)=-diff(interp1(r,Phi(1,:),[r1,r2],'spline'));
    rho(i)=pi*(L(i)^2-D^2)/D.*(DV(i)./FQ(1,1));
end
figure;
semilogx(L,rho);
