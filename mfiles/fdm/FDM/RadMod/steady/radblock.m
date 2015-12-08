function [Phi,Q,Psi]=radblock(r,z,kr,kz,FH,FQ)
% [Phi,Q,QRF,QFF]=radblock(r,z,kr,kz,FH,FQ)
% block centred, steady-state radial finite difference groundwater model
% r is radial node coordiantes (row vector, 1xNr+1)
% z is vertical coordinates (column vector, Nzx1+1)
% kr and kz radial and vertical conductivities (Nz)x(Nr)
% FH fixed heads (Nx*Nr), NaN for ordinary nodes, no NaN=fixed heads
% FQ nodal flows (Nz*Nr)
% Phi computed nodal heads
% Q   computed nodal flows
% Psi is srream function, the base of the aquifer is assumed impervious
% TO 990521 010811, 060405

%test data
if nargin==0,
	r=logspace(0,4,30);	Nr=length(r); rm=0.5*(r(1:end-1)+r(2:end));
	z=[0:-5:-50];		Nz=length(z); dz=abs(diff(z)); zm=0.5*(z(1:end-1)+z(2:end));
    ka=25; kc=0.01;
    kr=ka*ones(Nz-1,Nr-1); kr(1,:)=kc;
	FH=NaN*zeros(Nz-1,Nr-1); FH(1,:)=0; FH(2:3,1)=-1;
    FQ=zeros(size(FH));
    [Phi,Q,Psi]=radblockctrd(r,z,kr,kr,FH,FQ);
    figure
    contour(rm,zm,Phi,50,'b'); hold on
    contour(r ,z ,Psi,50,'r');
    QW=sum(Q(1:end,1));
    lambda=sqrt(ka*sum(dz(2:end))*(0.5*dz(1)/kc));
    fi=QW/(2*pi*ka*sum(dz(2:end)))*besselk(0,r/lambda);
    figure
    semilogx(rm,Phi,'r',r,fi,'b');
    return
end


HUGE=1e20;

FQ=FQ(:);
FH=FH(:);  Isfixhd=~isnan(FH); FH(find(isnan(FH)))=0;

r=r(:)'; dr=abs(diff(r)); Nr=length(r)-1; rm=0.5*(r(1:end-1)+r(2:end));
z=z(:) ; dz=abs(diff(z)); Nz=length(z)-1; zm=0.5*(z(1:end-1)+z(2:end));
Nod=Nr*Nz;

% node numbering
Nodes = reshape([1:Nod],Nz,Nr);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

r2=r(2:end); r1=r(1:end-1);

Rr=1/(2*pi)*(...
   1./kr(:,1:end-1).*((1./dz)*log(r( 2:end-1)./rm(1:end-1)))+...
   1./kr(:,2:end  ).*((1./dz)*log(rm(2:end  )./r( 2:end-1)))...
   );

Rz=0.5*dz*(1./(r(2:end).^2-r(1:end-1).^2))./kz/pi;
Rz=Rz(1:end-1,:)+Rz(2:end,:);
   
er = 1./Rr;
ez = 1./Rz;

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [er(:);er(:);ez(:);ez(:)],Nod,Nod,5*Nod);
Adiag=-sum(A,2);

Phi=spdiags(Adiag+Isfixhd*HUGE,0,A)\(HUGE*Isfixhd.*FH+~Isfixhd.*FQ);

Q=spdiags(Adiag,0,A)*Phi;

Phi=reshape(Phi,Nz,Nr);
Q  =reshape(Q,  Nz,Nr);
Psi=zeros(Nz+1,Nr+1); Psi(1:end-1,2:end-1)=flipud(cumsum(flipud(diff(Phi,1,2).*er)));
