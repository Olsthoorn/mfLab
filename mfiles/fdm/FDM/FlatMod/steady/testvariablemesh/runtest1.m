% for the modeflow keynote
% regular mesh, comparing meshctred with block centred and analytical
% TO 001019

% two layer model
% for the modeflow keynote in 2001
% comparing meshctrd with blockcntrd and analytical
% as in the thesis
% TO 001019
% TO 001225

close all
n=0.001;
L=600;

x=[0,logspace(0,log10(L),10)];
%x=[0:L/10:L];
y=[0,-1,-2,-3]';
% Block centerd model
x1=[-0.01,x,L+0.01];										% add small cells to ends to fix boundaries

K=[10;1e-3;30];

xm=0.5*(x1(1:end-1)+x1(2:end));							% cell centers
y1=[0,-1,-2,-3]';
dx1=diff(x1);													% cell width
dy1=abs(diff(y1));
Nx=length(dx1);												% number of cells in a row
Ny=length(dy1);												% number of cells in a column
k1=K.*dy1*ones(size(dx1));											% cell conductivities
FH1=NaN*ones(Ny,Nx); FH1([1,3],1)=0; FH1([1,3],end)=0;		% fixed heads
FQ1=zeros(size(FH1));
FQ1(1,:)=dx1*n;												% fixed flows (precipitation flow per cell)
[Phi1,Q1]=flatblockctrd(x1,y1,k1,k1,FH1,FQ1);		% solution for cells


% Mesh centerd grid
x2=x;
y2=y1;
dx2=diff(x2);
dy2=abs(diff(y2));
Nx2=length(x2);
Ny2=length(y2);
k2=K.*dy2*ones(size(dx2));									% element conductivities
FH2=NaN*ones(Ny2,Nx2); FH2(:,1)=0; FH2(:,end)=0;	% fixed heads
FQ2=zeros(size(FH2));
FQ2(1,:)=n*([dx2,0]+[0,dx2])/2;
[Phi2,Q2]=flatmeshctrd(x2,y2,k2,k2,FH2,FQ2);			% solution

%analytisch
ka=K(1);
xa=[0:L/50:L];
h=-n*((xa-L/2).^2-(L/2).^2)/(2*ka*dy1(1));


kD=K([1,3]).*dy1([1,3]);   kD=[kD,kD,kD,kD];
c =[1e6;dy1(2)/K(2)];            c=[c,c,c,c];
h=[0,0,0,0];
N=[0,n,n,0];
xa=[x(1),x(end)/2,x(end)];
X=[x(1):10:x(end)];

Q=[0,0,0;0,0,0];
[Section,Phi0,q,s]=NsecN(xa,kD,c,h,N,Q,X);

Q=[1,0,1;0,0,0];
[Section,Phia,q,s]=NsecN(xa,kD,c,h,N,Q,X);

Q=[0,0,0;1,0,1];
[Section,Phib,q,s]=NsecN(xa,kD,c,h,N,Q,X);

qq=[Phia(1,1)-Phi0(1,1),Phib(1,1)-Phi0(1,1);Phia(2,1)-Phi0(2,1),Phib(2,1)-Phi0(2,1)]\[-Phi0(1,1);-Phi0(2,1)];

Q=[qq(1),0,qq(1);qq(2),0,qq(2)];
[Section,Phi,q,s]=NsecN(xa,kD,c,h,N,Q,X);

figure;
plot(xm',Phi1([1,3],:)','r');	% Block centerd solution (leaving out the heads in the resistance layer)
hold on
plot(x2',Phi2','b');				% mesh centerd solution
%plot(xa',h'   ,'g');			% one layer analytic solution
plot(X',Phi'  ,'m');				% Analytic solution


figure;
plot(xm',Phi1([1,3],:)','+r',x2',Phi2([1,3],:)','xb',X',Phi','m');
legend('BC1','BC2','MC1','MC2','A1','A2');

hold on
ylim=get(gca,'ylim');
for i=1:length(x2)
   plot([x2(i),x2(i)],ylim);
end

xlabel('x [m]'); ylabel('head [m]');
title('comparison block- and mesh-centerd FD model');