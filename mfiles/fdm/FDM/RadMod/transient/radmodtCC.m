function [Phi,APhi,QR,QZ,St]=radmodtCC(r,z,t,kr,kz,ss,FH,FQ,Phi0,theta)
% [Phi,APhi,QR,QZ,St]=RadModCC(r,z,t,kr,kz,ss,FH,FQ,Phi0,theta)
% Cell-centred, transient radial finite difference groundwater model
% r is radial NCelle coordiantes (row vector, 1xNr)
% z is vertical coordinates (column vector, Nzx1)
% t times to compute results
% kr and kz radial and vertical conductivities (Nz)x(Nr)
% ss storativity [1/m] (Nz)*(Nr)
% FH fixed heads (Nx*Nr), NaN for ordinary Cells, no NaN=fixed heads
% FQ NCellal flows (Nz*Nr)
% Phi0 = matrix with initial Phi
% theta degree of implicitness, use 0.5 if possible or 0.67 must be between 0 and 1
% Phi computed Cell heads
% APhi= computed Cell flows Q, the water leaving node through the model computed as A*Phi
% The water balance is APhi=FQ-St to get nodal flows compute APhi+St this includes fixed head nodes
% St  = Storage in cells during the time step
% QR  = Q right face during time step
% QZ  = Q upward z-direction across face during time step
% TO 990521 051221

HUGE=1e20;

FH=FH(:); FQ=FQ(:); Phi0=Phi0(:);
Isfixhd=~isnan(FH(:)); FH(find(isnan(FH)))=0;

r=r(:)'; dr=abs(diff(r)); Nr=length(dr);  rm=0.5*(r(1:end-1)+r(2:end));
z=z(:) ; dz=abs(diff(z)); Nz=length(dz);
t=t(:);  dt=abs(diff(t)); Nt=length(t);
NCell=Nr*Nz;

R=ones(size(dz))*r; Rm=ones(size(dz))*rm;

% Cell numbering
Cells = reshape([1:NCell],Nz,Nr);
Il=Cells(:,2:end);   Jl=Cells(:,1:end-1);
Ir=Cells(:,1:end-1); Jr=Cells(:,2:end);
It=Cells(2:end,:);   Jt=Cells(1:end-1,:);
Ib=Cells(1:end-1,:); Jb=Cells(2:end,:);
Il=Il(:); Jl=Jl(:); Ir=Ir(:); Jr=Jr(:); It=It(:); Jt=Jt(:); Ib=Ib(:); Jb=Jb(:);

%Matrix coefficients
Tr=(dz*ones(size(dr))).*kr;
re=(log([Rm(:,2:end-1),R(:,end)]./R(:,2:end-1))./Tr(:,2:end)+...            % includes whole of last column
    log(R(:,2:end-1)./[R(:,1),Rm(:,2:end-1)])./Tr(:,1:end-1))./(2*pi);      % includes whole of first column
er=1./re;

rz=(dz*ones(1,Nr))./(kz.*pi.*(R(:,2:end).^2-R(:,1:end-1).^2))/2;
rz([1,end],:)=rz([1,end],:)*2.0;                                            % to include whole of top and bottom layer
rz=0.5*(rz(2:end,:)+rz(1:end-1,:)); 
ez=1./rz;
 
es= pi*ss.*(dz*(r(2:end).^2-r(1:end-1).^2));
 
%Maxtrix construction
A=-sparse([Il;Ir;It;Ib],[Jl;Jr;Jt;Jb],[er(:);er(:);ez(:);ez(:)],NCell,NCell,5*NCell);
Adiag=-sum(A,2);

% initiate matrices
Phi(1:NCell,1)=Isfixhd.*FH+~Isfixhd.*Phi0;  % store initial heads
dt=diff(t);
QR(1:Nz,1:Nr+1,1)=0;
QZ(1:Nz+1,1:Nr,1)=0;
St(1:Nz*Nr,1)=0;
APhi(1:Nz*Nr,1)=0;

% loop over time steps
for it=2:Nt;
	ES=es(:)/(dt(it-1)*theta);
	PhiM=spdiags(Adiag+Isfixhd*HUGE+~Isfixhd.*ES,0,A)\...
   	        (HUGE*Isfixhd.*FH+~Isfixhd.*(FQ+ES.*Phi(:,it-1)));
    Phi (:,it  )= Phi(:,it-1)+(PhiM-Phi(:,it-1))/theta;
    St  (:,it-1)=(Phi(:,it)-Phi(:,it-1)).*es(:);
	APhi(:,it-1)= spdiags(Adiag,0,A)*PhiM;

    PhiM=reshape(PhiM,Nz,Nr);
    QR(:,2:end-1,it-1)=er.*(PhiM(:,1:end-1)-PhiM(:,2:end));
    QZ(2:end-1,:,it-1)=ez.*(PhiM(2:end,:)-PhiM(1:end-1,:));
end
Phi =reshape(Phi,Nz,Nr,Nt);    % head per cell m
APhi=reshape(APhi,Nz,Nr,Nt-1);    % flow into cell (source) [L3/T]
St  =reshape(St  ,Nz,Nr,Nt-1);