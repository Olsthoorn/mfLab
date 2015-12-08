% radblockctrd test 2, vlakke stroming met bouwput en onttrekking
% berekende met radiaal blok gecentreerde model dat ook stroomfunctie berekend
% (c) Olsthoorn, 000523
% TO 000530, stroomfunctie ingebouwd. Is nu in orde.

% Opzetten van het radiale model van een bouwput met een damwand waar, tegen de binnenwand een onttrekkings
% scherm wordt geplaatst met een vast stijghoogte PhiW. Het model is semispanningswater. De stijghoogte boven
% de slecht doorlatende laag is Phi0. kDW is de k van de DamWand in resp. r- en z-richting.
c=1000; kr=30; kz=30; kDWr=1e-3; kDWz=1e-3;
PhiW =-3; Phi0 = 0;

RBouwput =18;		Ddamwand=0.10;		R=1e4;

% Afmetingen: hoogten van de kenmerkende vlakken van het model 
z0   =0;				% bovenzijde toplaag
z1   =-10;			% onderzijde toplaag
zWtop=-11;			% top onttrekkingsscherm
zWbot=-16;			% onderkant onttrekkingsscherm
zDW  =-18;			% onderkant damwand
zBot =-70;			% onderkant WVP

ii=10;

% Even een geschikt aantal r-waarden aanmaken voor het modelgrid (zie help logspace)
% we zorgen voor een fijnmazig netwerk in de buurt van de damwand.
r=[RBouwput-fliplr(logspace(-1,log10(RBouwput),20)),RBouwput,RBouwput+logspace(log10(Ddamwand),log10(R),20)];

% een geschikt aantal z-waarden aanmaken voor het modelgrid. We zorgen voor een heel fijnmazig netwerk
% ter hoogte van de tip van de damwand
z=[z0;z1;[zWtop:-0.25:zWbot]';[zWbot-0.25:-0.25:zDW]';[zDW-logspace(-1,1,10)]';[zDW-12:-5:zBot]'];

r(1)=0.1;		% minimum radius van het model, voorkomt delen door nul

Nr=length(r); Nz=length(z);			% Tellen van aantal cellen in r- en z-richting

rm=(r(1:end-1)+r(2:end))/2;			% radii van middens van de cellen
zm=(z(1:end-1)+z(2:end))/2;			% hoogte van de middens van de cellen
FH=NaN*zeros(Nz-1,Nr-1);				% declareren (aanmaken) matrix voor stijghoogterandvoorwaarden
FQ=    zeros(Nz-1,Nr-1);				% idem voor gegeven volumestromen

Kr=kr*ones(Nz-1,Nr-1);					% matrix horizontale doorlatendheden
Kz=kz*ones(Nz-1,Nr-1);					% matrix verticale doorlatendheden
Kr(1,:)=abs(z(2)-z(1))/(2*c);			% verander k toplaag in waarde die de weerstand creert
Kz(1,:)=abs(z(2)-z(1))/(2*c);			% idem verticale doorlatendheid

I=min(find(rm>RBouwput));				% zoek in welke cellen de damwand zit
J=find(zm>zDW);							% idem in verticale richting
for i=1:length(I)							% voor elke kolom met damwand
  Kr(J,I(i))=kDWr;						% vervang kr door die van damwand
  Kz(J,I(i))=kDWz;						% vervang kz door die van damwand
end

I=max(find(r<RBouwput));				% we gaan hier de het onttrekkingsscherm tegen de binnenzijde damwand plaatsen
J=find(zm>=zWbot & zm<=zWtop);		% zoek cellen waar de stijghoogte van het scherm moet worden opgelegd
for i=1:length(I)							% voor alle kolommen met een scherm
   FH(J,I(i))=PhiW;						% stop deze stijghoogte in de randvoorwaarde matrix voor alle rijen
end
Iw=find(~isnan(FH));						%
FH(1,:)=0;

[Phi,Q,Psiy,Psix]=flatblockctrd(r,z,Kr,Kz,FH,  FQ);
%[Phi,Q,Psiy,Psix]=radblockctrd( r,z,Kr,Kz,FH,  FQ);			% run het model
% Het resultaat is de stijghoogte in de cellen de berekende infiltratie in elke cel, de stroomfunctie
% van de celwanden, met verticale cuts en idem met horizontale cuts

Qw=sum(Q(Iw))								% totale onttrekking (infiltratie positief)

% bereken een geschikte range voor het tekenen van contourlijnen
m=min([Psiy(:);Psix(:)]); M=max([Psiy(:);Psix(:)]); d=(M-m)/30; psirange=[m:d:M];		% stroomfuncie
m=min(Phi(:)); M=max(Phi(:)); d=(M-m)/120; phirange=[m:d:M];										% stijghoogten

figure;												% contouren stijghoogten en stroomfunctie met horizontale cuts
contour(rm,zm,Phi,phirange); hold on		% contouring stijghoogte
contour(r,z,Psix,psirange);					% contouring met horizontale cuts
% let op: De stroomfunctie is alleen bekend op de knooppunten, de stijghoogte in de cellen

xlabel('r [m]'); ylabel('z [m]'); title('flatblockctrd, Phi en Psix');
set(gca,'ylim',[zBot,0],'xlim',[0,abs(zBot)]);


figure;												% contouren stijghoogten en stroomfunctie met verticale cuts
contour(rm,zm,Phi,phirange); hold on		% contouring stijghoogte
contour(r,z,Psiy,psirange);					% contouring stroomfuncties met verticale cuts
% let op: De stroomfunctie is alleen bekend op de knooppunten, de stijghoogte in de cellen

xlabel('r [m]'); ylabel('z [m]'); title('flatblockctrd, Phi en Psiy');
set(gca,'ylim',[zBot,0],'xlim',[0,abs(zBot)]);

% analytische verlaging volgens de Glee met bouwputonttrekking ter vergelijking met numerieke model
kD=sum(Kr(2:end,1).*abs(diff(z(2:end))));		% Even de kD uitrekenen uit de k-bijdragen van de modellagen
phi =Qw/(2*pi*kD)*besselk(0,r/sqrt(kD*c));	% De Glee

%figure;												% Tekenen verlaging volgens de Glee en model (boven en ondervlak WVP)
%plot(r,phi,'r'); hold on;
%plot(rm,Phi([2,end] ,:),'b');
%set(gca,'xscale','log');
%xlabel('r [m]'); ylabel('phi [m]'); title('analytisch (De Glee) vs. radmod');