function [vr,vz]=radmodtCCvelocity(r,z,t,kr,kz,por,QR,QZ)
% [vrvz]=radmodtCCvelocity(r,z,t,por,QR,QZ)
% Computes velocity field (quiverfield)
% r is radial NCelle coordiantes (row vector, 1xNr)
% z is vertical coordinates (column vector, Nzx1)
% t times to compute results
% kr and kz radial and vertical conductivities (Nz)x(Nr)
% por porosity (Nz)*(Nr)
% QR  = Q right face during time step
% QZ  = Q upward z-direction across face during time step
% TO 051222

r=r(:)'; dr=abs(diff(r)); Nr=length(dr);  rNode=0.5*(r(1:end-1)+r(2:end)); rNode(1,end)=r(1,end);
z=z(:) ; dz=abs(diff(z)); Nz=length(dz);  zNode=0.5*(z(1:end-1)+z(2:end)); zNode(1,end)=z(1,end);
t=t(:);  dt=abs(diff(t)); Nt=length(t); Ndt=length(dt);

ARNode=2*pi.*por.*(dz*rNode);        % total cross section at cell (2 pi r Dz)
AZNode=  pi.*por.*(ones(size(dz))*(r(:,2:end).^2-r(:,1:end-1).^2));      % total horizontal cross section

vr=zeros(Nz,Nr,Ndt);
vz=zeros(Nz,Nr,Ndt);

for idt=1:Ndt
    vr(:,:,idt)=0.5*(QR(:,1:end-1,idt)+QR(:,2:end,idt))./ARNode;
    vz(:,:,idt)=0.5*(QZ(1:end-1,:,idt)+QZ(2:end,:,idt))./AZNode;
end
