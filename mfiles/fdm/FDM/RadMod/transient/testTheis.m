% radmodt test
% vergelijking radmodt en Theis
% TO 990522

r=logspace(-1,2,31);
y=[0,-100]';
D=y(1)-y(2);

Nr=length(r);
Ny=length(y);

kr=10;
ky=10;
ss=1e-3;

t=logspace(-3,-1,41)';
u=ss*D/(4*kr*D)*r.^2/t(1);
FQ  =zeros(Ny,Nr); FQ(:,1)=-1200;
phiTh=sum(FQ(:,1))/(4*pi*kr*D)*expint(u);

Phi0=[phiTh;phiTh];											% initialize with theis(r,t(1));
FH  =zeros(Ny,Nr)*NaN;										% no fixed heads
FQ  =zeros(Ny,Nr); FQ(:,1)=-1200;						% fixed extractions

theta=0.5;
[Phi,Q]=radmodt(r,y,t,kr,ky,ss,FH,FQ,Phi0,theta);

u=ss*D/(4*kr*D)*(ones(size(t))*r.^2)./(t*ones(size(r)));
phiTh=sum(FQ(:,1))/(4*pi*kr*D)*expint(u);

figure; semilogx(r,squeeze(Phi(end,:,:)),r,phiTh(:,:),'+');
xlabel('r [m]'); ylabel('s [m]'); title('vergelijking model (-) en Theis (+)');

figure; semilogx(t,squeeze(Phi(1,:,:)),t,phiTh(:,:),'+');
xlabel('t [d]'); ylabel('s [m]'); title('vergelijking model (-) en Theis (+)');