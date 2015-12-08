function Out=wideDitch1t(tne,P)
%Out=wideDitch1t(tne,P)
%
%
error('Dont''use, not yet tested !!');
%
% Transient wideditch head and seepage computation
%
%
% Analytic solution for flow in shallow aquifer on top of semi-confined
% regional aquifer with given head PHI.
% Situation
%
%                               |||||||||||| N2 |||||||||||||||||||
%     wide ditch ,h=hLR      |  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv |
%                            |                                      |
%    |~~~~~~~~  csb~~~~~~~~~~~                                      |
%    |          kD1                         kD2                     |
%    |//////////c1//////////////////////////c2///////////////////// |
%    |          PHI                         PHI                     |
%    |^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ |
%    |||||||||||||q1|||||||||||||||||||||||  q2 ||||||||||||||||||| |
%    |==== x1 ===>                                   <== x2 ========|
%    |<-------- b1 --------->|<-------------- b2 ------------------>|
%
% Ditch can be on either side, like the recharge.
%
% Flux through confining bed positive if upward (Dutch "kwel").
%
%
% TO 101227 110122

%% Initialize solution specific
Np=length(P);

t=tne(:,1);
Nt = length(t);
Prec=tne(:,2)-tne(:,3);
Dt=diff(t); Dt=[Dt(1); Dt];

%% 

dd = 0.3;  % [ m ] ditch depth
dw = 2;    % [ m ]ditch width
fb = 0.6;  % [ m ] freeboard (drooglegging)
csb = ones(size(P))*[ 1  Inf];  % [ m ] ditch bottom resistance, use Inf if no ditch

PHI = [P.phi]';               % [ m ] prescribed head in regional aquifer
hLR = [P.AHN]'-0.6;        % [ m ] head in ditch
kD  = [ [P.hk1]'.*[P.D1]'-fb-dd   [P.hk1]'.*[P.D1]'-fb];

C   = [ [P.c]' [P.c]']; % [ d ] resistance at bottom of cover layer
b   = [ dw/2*ones(size(P)) [P.b]'-dw/2]; % [ m ] with of the two compartments

c   = 1./(1./csb+1./C); % [ d ] net resistance in compartments
lam = sqrt(kD.*c);      % [ m ] characteristic length of the compartments
mu  = [P.mu]';

B=NaN(length(P),2);             % [ - ] dimensionless constants
B(:,1)=cosh(b(1)./lam(1))+(kD(1)./lam(1).*sinh(b(1)./lam(1)))./(kD(2)./lam(2).*sinh(b(2)./lam(2))).*cosh(b(2)./lam(2));
B(:,2)=cosh(b(2)./lam(2))+(kD(2)./lam(2).*sinh(b(2)./lam(2)))./(kD(1)./lam(1).*sinh(b(1)./lam(1))).*cosh(b(1)./lam(1));


phi = (hLR./csb+PHI./C)./(1./csb+1./C);

%% Transient

beta2=B(:,2)./(lam(:,2)./b(:,2).*sinh(b(:,2)./lam(:,2)));
T=mu.*c(:,2);

hm=NaN(Np,Nt);

e=exp(-Dt(1)./T);
theta=Prec(1)*c(:,2)+(Prec(1).*c(:,2)+phi(:,2)-phi(:,1))./beta2 +phi(:,2)./c(:,2);
hm(1)=phi(:,2)+theta*(1-e);

for it=2:Nt
    e=exp(-Dt(it)./T);
    theta=Prec(it)*c(:,2)+(Prec(it)*c(:,2)+phi(:,2)-phi(:,1))./beta2 +phi(:,2)/c(:,2);
    hm(:,it)=phi(:,2)+(hm(:,it-1)-phi(:,2))*e+theta*(1-e);
end

%% Plot results
iP=1;

figure; hold on; grid on; xlabel('x [m]'); ylabel('head [m]');
title(sprintf('wide ditch, head versus time section %d q=%.3f m/d',iP,q));
leg='';

plot(t,hm(iP,:),'r');

legend(leg)

Out=hm;
