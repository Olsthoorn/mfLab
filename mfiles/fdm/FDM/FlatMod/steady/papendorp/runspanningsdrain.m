%rundrain
%berekent Q door drain in semi-spanningswater met gemiddelde verlaging s0
%berekent bovendien q en qinf voor oneindig lange drain en de verhouding tussen beide
%
% TO 001006
close all

s0=2.1;			% doelverlaging

% gegevens pakket en situatie, drain in oneindig semi-spanningspakket
L=40;				% lengte drain
n=100;			% aantal stapjes langs drain
D=50;				% Pakketdikte
k=45;				% doolatendheid
kD=k*D;			% kD
c=500;			% c
Lambda=sqrt(kD*c); % spreidingslengte
B=3;				% breedte drain

% Eerst berekenen verlaging voor een dummy Q van 1 m3/d
Q=1;				% dummy startwaarde
q=Q/L;			% specifiek, nog dummy
v=q/B;			% snelehid,  nog dummy
[x,s]=semispanningsdrain(q,L,n,kD,c);	% integreert bijdrage van drainstukjes langs de hele drain gegeven q
plot(x,s);
title('dummy resultaat, voor Q=1');

% Berekening van de Q bij de doelverlaging, opschaling Q
sM=mean(s,2);		% gemiddelde s volgens dummy
Q=s0/sM*Q;			% opschalen naar doel s, s0
q=Q/L;				% nieuwe q
v=q/B;				% nieuwe v
s=s0/sM*s;			% opschalen naar doel s
qinf=2*s0*kD/Lambda;	% q voor oneindig lange sleuf, zonder intredeweerstand
figure
plot(x,s);			% plotten situatie bij gemiddelde verlaging s0
title(sprintf('Resultaat,s=%.2f, Q=%.1f m3/d, q=%.2f m2/d, qinf=%.2f m2/d n=q/qinf=%.1f',s0,Q,q,qinf,q/qinf));

% partial penetration
%weerstand pp zou zijn:
cpp=B/(pi*k)*log(D/B);		% weerstand pp
dspp=cpp*Q/L/B;				% verlaging door deze weerstand in geval Q=Q zoals hiervoor berekend
frac=s0/(dspp+s0);			% verlaging terugschalen totdat totale verlaging gelijk is aan s0, dit is de schaal

s=frac*s;						% verlagingsverloop terugschalen
Q=frac*Q;						% onttrekking terugschalen
q=Q/L; v=q/B; dspp=cpp*v;	% idem q en v
figure							% laat maar zien
plot(x,s);						% laat maar zien
title(sprintf('Resultaat,s=%.2f,dspp=%.2f,Q=%.1f m3/d,q=%.2f m2/d,qinf=%.2f m2/d,n=q/qinf=%.1f',mean(s),dspp,Q,q,qinf,q/qinf));
