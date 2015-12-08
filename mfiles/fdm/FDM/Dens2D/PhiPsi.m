function [Phi,Q,PsiStruct,PsiJump,Psi]=phipsi(x,z,kx,kz,FH,FQ,NPhi,NPsi);
% PHIPSI(x,z,kx,kz,FH,FQ,NPsi,NPsi);
% computes and plots the head Phi and stream function Psi for plane groundwater flow
% using mesh-centred finite differences.
% x,z are x and z node coordinate vectos resp.
% kx and ky the conductivities
% FH the fixed heads Nx*Ny matrix of NaN's and fixed heads at given nodes
% FH is Nz*Nx matrix of nodal injections (extractions negative)
% NPhi is number of head contours, NPsi is number of Psi contours.
% PsiStruct is the Psi but plit where multivalued. use contour(PsiStruct(i).x,z,PsiStruct(i).Psi); to plot.
% Uses flatmeshctrd and psimeshctrd which includes the internal sources according to Van den Akker (1982).
% It assumes vertical branchcuts.
% TO 000528

[Phi,Q]=flatmeshctrd(x,z,kx,kz,FH,FQ);

dddy=zeros(size(FH));		% density may be included in future

[Psi,Rot,FPsi,PsiJump]=psimeshctrd(x,z,kx,kz,Q,dddy);

[phirange,dphi]=nicerange(Phi,NPhi);
[psirange,dpsi]=nicerange(Psi,NPsi);

% find internal sources:
split=any(FQ(2:end,:)~=0|~isnan(FH(2:end,:))); split([1,end])=[0,1]; I=find(split);
i1=1;
for i=1:length(I)
   i2=I(i);
   PsiStruct(i).Psi=Psi(:,i1:i2);
   PsiStruct(i).Psi(:,1)=PsiStruct(i).Psi(:,1)-PsiJump(:,i1);
   PsiStruct(i).x=x(i1:i2);
   i1=i2;
end

figure
contour(x,z,Phi,phirange,'b');
hold on
for i=1:length(PsiStruct)
   contour(PsiStruct(i).x,z,PsiStruct(i).Psi,psirange,'r');
end
s=sprintf('PhiPsi, dPhi=%%.%df m, dPsi=%%.%df m2/d',abs(min(0,fix(log10(dphi))-1)),abs(min(0,fix(log10(dpsi))-1)));
xlabel('x [m]'); ylabel('z [m]'); title(sprintf(s,dphi,dpsi));
