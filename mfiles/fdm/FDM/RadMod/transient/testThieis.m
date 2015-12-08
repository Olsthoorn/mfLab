% radmodtest
% thiem
r=logspace(-1,3,20);
y=[0,-1]';

Nr=length(r);
Ny=length(y);

kr=10;
ky=10;
ss=1e-3;

Phi0=zeros(Ny,Nr);
FH  =zeros(Ny,Nr)*NaN;
FQ  =zeros(Ny,Nr);
FH(:,1)=-1;
t=logspace(-3,3,61);

[Phi,Q]=radmod(r,y,t,kr,ky,ss,FH,FQ,Phi0,theta);
semilogx(r,Phi(end,:,:)
