function [Phi,Q]=radmodt(r,z,t,kr,kz,ss,FH,FQ,Phi0,theta)
% [Phi,Q]=RadMod(r,z,t,kr,kz,ss,FH,FQ,Phi9,theta)
% mesh centred, transient radial finite difference groundwater model
% r is radial node coordiantes (row vector, 1xNr)
% z is vertical coordinates (column vector, Nzx1)
% t times to compute results
% kr and kz radial and vertical conductivities (Nz-1)x(Nr-1)
% ss storativity [1/m] (Nz-1)*(Nr-1)
% FH fixed heads (Nx*Nr), NaN for ordinary nodes, no NaN=fixed heads
% FQ nodal flows (Nz*Nr)
% Phi0 = matrix with initial Phi
% theta degree of implicitness, use 0.5 if possible or 0.67 must be between 0 and 1
% Phi computed nodal heads
% Q   computed nodal flows
%
% TO 990521

HUGE=1e20;

FH=FH(:); FQ=FQ(:); Phi0=Phi0(:);
Isfixhd=~isnan(FH(:)); FH(find(isnan(FH)))=0;

r=r(:)'; dr=abs(diff(r)); Nr=length(r);
z=z(:) ; dz=abs(diff(z)); Nz=length(z);
t=t(:);  dt=abs(diff(t)); Nt=length(t);
Nod=Nr*Nz;

Nkr=size(kr,1);
if size(kr,2)<Nr-1,   kr=[kr,kr(:,end)*ones(1,Nr-1-size(kr,2))]; end
if size(kr,1)<Nz-1,   kr=[kr;ones(Nz-1-size(kr,1),1)*kr(end,:)]; end
if size(kz,2)<Nr-1,   kz=[kz,kz(:,end)*ones(1,Nr-1-size(kz,2))]; end
if size(kz,1)<Nz-1,   kz=[kz;ones(Nz-1-size(kz,1),1)*kz(end,:)]; end
if size(ss,2)<Nr-1,   ss=[ss,ss(:,end)*ones(1,Nr-1-size(ss,2))]; end
if size(ss,1)<Nz-1,   ss=[ss;ones(Nz-1-size(ss,1),1)*ss(end,:)]; end

% node numbering
Nodes = reshape([1:Nod],Nz,Nr);
Il=Nodes(:,2:end);   Jl=Nodes(:,1:end-1);
Ir=Nodes(:,1:end-1); Jr=Nodes(:,2:end);
It=Nodes(2:end,:);   Jt=Nodes(1:end-1,:);
Ib=Nodes(1:end-1,:); Jb=Nodes(2:end,:);
Il=Il(:); Jl=Jl(:); Ir=Ir(:); Jr=Jr(:); It=It(:); Jt=Jt(:); Ib=Ib(:); Jb=Jb(:);

r2=r(2:end); r1=r(1:end-1);

er = 2*pi*kr.*(dz*ones(size(dr)))./log(ones(size(dz))*(r2./r1));
er = 0.5*([zeros(size(dr));er]+[er;zeros(size(dr))]);

rMean=(r2-r1)./log(r2./r1);                % gradient at rMean complies with ordinary Darcy
alfa1=(rMean.^2-r1.^2)./(r2.^2-r1.^2);
alfa2=(r2.^2-rMean.^2)./(r2.^2-r1.^2);

ez= kz.*pi.*(ones(size(dz))*(r2.^2-r1.^2))./(dz*ones(size(dr)));
ez= (ones(size(dz))*[0,alfa2]).*[zeros(size(dz)),ez]+...
    (ones(size(dz))*[alfa1,0]).*[ez,zeros(size(dz))];
 
es= pi*ss.*(dz*(r2.^2-r1.^2));
es= (ones(size(dz))*[0,alfa2]).*[zeros(size(dz)),es]+...
    (ones(size(dz))*[alfa1,0]).*[es,zeros(size(dz))];
es=0.5*([zeros(size(r));es]+[es;zeros(size(r))]);
 
er=er(:); ez=ez(:); es=es(:);
 
A=-sparse([Il;Ir;It;Ib],[Jl;Jr;Jt;Jb],[er;er;ez;ez],Nod,Nod,5*Nod);
Adiag=-sum(A,2);

Phi(1:Nod,1)=Isfixhd.*FH+~Isfixhd.*Phi0;
  Q(1:Nod,1)=FQ;
dt=diff(t);
for it=2:Nt;
	ES=es/(dt(it-1)*theta);
	PhiM=spdiags(Adiag+Isfixhd*HUGE+~Isfixhd.*ES,0,A)\...
   	        (HUGE*Isfixhd.*FH+~Isfixhd.*(FQ+ES.*Phi(:,it-1)));
%   Phi(:,it)=Isfixhd.*FH+~Isfixhd.*(Phi(:,it-1)+(PhiM-Phi(:,it-1))/theta);
    Phi(:,it)=Phi(:,it-1)+(PhiM-Phi(:,it-1))/theta;

	QM=spdiags(Adiag,0,A)*PhiM;
	Q(:,it)=QM;
end
Phi=reshape(Phi,Nz,Nr,Nt);
Q  =reshape(Q  ,Nz,Nr,Nt);
