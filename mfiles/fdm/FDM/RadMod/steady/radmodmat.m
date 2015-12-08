function [Phi,Q]=radmod(r,z,kr,kz,FH,FQ)
% [Phi,Q]=RadMod(r,z,kr,kz,FH,FQ)
% mesh centred, steady-state radial finite difference groundwater model
% r is radial node coordiantes (row vector, 1xNr)
% z is vertical coordinates (column vector, Nzx1)
% kr and kz radial and vertical conductivities (Nz-1)x(Nr-1)
% FH fixed heads (Nx*Nr), NaN for ordinary nodes, no NaN=fixed heads
% FQ nodal flows (Nz*Nr)
% Phi computed nodal heads
% Q   computed nodal flows
%
% TO 990521

HUGE=1e20;

FQ=FQ(:);
FH=FH(:);  Isfixhd=~isnan(FH); FH(find(isnan(FH)))=0;

r=r(:)';
z=z(:);

Nr=length(r);
Nz=length(z);
Nod=Nr*Nz;

Nkr=size(kr,1);
if size(kr,2)<Nr-1,   kr=[kr,kr(:,end)*ones(1,Nr-1-size(kr,2))]; end
if size(kr,1)<Nz-1,   kr=[kr;ones(Nz-1-size(kr,1),1)*kr(end,:)]; end
if size(kz,2)<Nr-1,   kz=[kz,kz(:,end)*ones(1,Nr-1-size(kz,2))]; end
if size(kz,1)<Nz-1,   kz=[kz;ones(Nz-1-size(kz,1),1)*kz(end,:)]; end

% node numbering
Nodes = reshape([1:Nr*Nz],Nz,Nr);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);

dr=ones(Nz-1,1)*diff(r);
dz=diff(z)*ones(1,Nr-1);

er0= pi*kr./log(r(1:end-1,2:end)./r(1:end-1,1:end-1)).*dz;
er= [er0;zeros(1,Nr-1)]+[zeros(1,Nr-1);er0];

ez1= kz.*pi.*((r(1:end-1,1:end-1)+dr/2).^2-(r(1:end-1,1:end-1)   ).^2)./dz;
ez2= kz.*pi.*((r(1:end-1,2:end)       ).^2-(r(1:end-1,2:end)-dr/2).^2)./dz;
ez = [ez1,zeros(Nz-1,1)]+[zeros(Nz-1,1),ez2];

A=-sparse([Il(:);Ir(:);It(:);Ib(:)],...
          [Jl(:);Jr(:);Jt(:);Jb(:)],...
          [er(:);er(:);ez(:);ez(:)],Nod,Nod,5*Nod);
Adiag=-sum(A,2);

Phi=spdiags(Adiag+Isfixhd*HUGE,0,A)\(HUGE*Isfixhd.*FH+~Isfixhd.*FQ);

Q=spdiags(Adiag,0,A)*Phi;
Phi=reshape(Phi,Nz,Nr);
Q  =reshape(Q,  Nz,Nr);