% run infiltratiegebieden
% This m file is a trial to run in much detail the groundwater model for the 2nd and 3rd recharge areas
% It now works, it only missis the drains. Coordinates are needed for that, plus entry resistance. But this
% is not a major problem
% It further misses inactive cells. These have not been built into Flatmod. They could be. Alternatively, just
% for this situation fix the heads outside the area, to no flow occurs underneath the canals, which is possible
% because we use and entry resistance (general head boundary) to represent the resistance of the surface water
%  to the aquifer but also for the connection with the lower aquifer. This connection has to be built in by the
% boundary condition, but is not a problem, because nothing needs to be changed to this model.
% Transient flow would be the last step. If small steps are taken, it is possible to compute the transient flow
% by an explicit Eulerian scheme, which is very fast, because no matrix needs to be solved. Only a transient
% initial situation would be benificial. An arbitrary initial situation can be used to see if the explicit scheme
% is stable and rapidly converges towards its steady state solution, which is likely with small steps.
% TO 010825
readdxf;			% reads dxf file of 2nd and third recharge area of AWD and loads characteristics of ponds and canals

clear
load PlineIGeb2&3

nn=0.003;

xM=0.5*(x(1:end-1)+x(2:end)); yM=0.5*(y(1:end-1)+y(2:end));
dx=diff(x); dy=diff(y); [DX,DY]=meshgrid(dx,dy);
c=cat(1,Owater{:,3});
A=DX.*DY;
N=A.*nn;
I=find(In~=0);

Entr=NaN*zeros(size(In));
Entr(I)=A(I)./c(In(I));
FH=0*Entr;
FQ=N;
kD=100*ones(size(In));
tic
[Phi,Q]=flatblockctrd(x,y,kD,kD,FH,FQ,Entr);
toc

figure
contour(xM,yM,Phi,[0:0.02:0.5]); xlabel('x'); ylabel('y'); title('Infiltratiegebieden II en III');
hold on;
for i=1:length(pline)
   plot(pline(i).x,pline(i).y,'b');
end
