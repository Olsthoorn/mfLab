% Voorbeeld van toepassing Radmod, De Glee, TO 990522
r=logspace(-5,4,91);     									% r op logschaal tussen 1 en 10000, 41 punten
D1=1; D2=50; yTop=1; y=yTop-[0,D1,D1+D2];				% y-waarden, afh van dikten
c=250;															% weerstand topsysteeem
kr=[D1/c;20]; ky=kr;											% weerstand toplaag en k laag 2
FH= zeros(length(y),length(r))* NaN; FH(1,:)=0;		% polderpeil
FQ=zeros(length(y),length(r)); FQ(2:end,1)= -1200;	% onttrekking, 2400 m3/d
[Phi,Q]=radmod(r,y,kr,ky,FH,FQ);							% run het model
Q0=sum(Q(1,:))													% check berekende infiltratie
contour(r,y,Phi);					% contour verkregen stijghoogten (hier niet zinvol)
PhiDeGlee=-Q0/(2*pi*kr(2,1)*D2)*besselk(0,r./sqrt(kr(2,1)*D2*c));  % DeGlee
figure; semilogx(r,Phi(end,:),r,PhiDeGlee,'r+');	grid; 	% grafiek, vergelijking met DeGlee
legend('model','analytisch, De Glee',4);
xlabel('r [m]'); ylabel('s [m]'); title('Test: Radmod vs De Glee')
