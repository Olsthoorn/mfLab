
solution_name='wideDitch';

%mktestset

% Steady-state wideditch head and seepage computation
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
    Q     =    0.3;         % [m2/d] extraction between the compartments
    qsoll = 0.01;         % [m/d] desired upward flux over entire width b1+b2
    N   =[ 0.004 0.004];        % [m/d] reharge on the two compartments
    PHI = -1;               % [ m ] prescribed head in regional aquifer
    hLR = [-0.5 -0.5];        % [ m ] head in ditch
    kD  = [  5  15];        % [m2/d] transmissivity of cover layer
    C   = [ 20  20];        % [ d ] resistance at bottom of cover layer
    b   = [ 5   45];        % [ m ] with of the two compartments
    csb = [ 1  1000000];        % [ m ] ditch bottom resistance, use Inf if no ditch
    w=1e6;      % [d/m] resistance through vertical side of ditch

hditch=hLR(1);  % [ m ] ditch level

c   = 1./(1./csb+1./C); % [ d ] net resistance in compartments
lam = sqrt(kD.*c);      % [ m ] characteristic length of the compartments

B=NaN(1,2);             % [ - ] dimensionless constants
B(1)=cosh(b(1)/lam(1))+(kD(1)/lam(1)*sinh(b(1)/lam(1)))/(kD(2)/lam(2)*sinh(b(2)/lam(2)))*cosh(b(2)/lam(2));
B(2)=cosh(b(2)/lam(2))+(kD(2)/lam(2)*sinh(b(2)/lam(2)))/(kD(1)/lam(1)*sinh(b(1)/lam(1)))*cosh(b(1)/lam(1));

beta=NaN(1,2);
beta(1)=lam(1)/b(1)*sinh(b(1)/lam(1))/B(1);
beta(2)=lam(2)/b(2)*sinh(b(2)/lam(2))/B(2);

dphidphi=1./(C./csb+1); % [ - ] dphi(1)/dphi and dphi(2)/dphi

dQdphi=NaN(1,2);
dQdphi(1)=(dphidphi(1)+(dphidphi(2)-dphidphi(1))*cosh(b(1)/lam(1))/B(1))/...
    (w+cosh(b(1)/lam(1))/B(1)/(kD(2)/lam(2)*tanh(b(2)/lam(2))));
dQdphi(2)=(dphidphi(2)+(dphidphi(1)-dphidphi(2))*cosh(b(2)/lam(2))/B(2))/...
    (w+cosh(b(2)/lam(2))/B(2)/(kD(1)/lam(1)*tanh(b(1)/lam(1))));

dhdphi=NaN(1,2);        % [ - ] dh(1)/dphi and dh(2)/dphi
dhdphi(1)=dphidphi(1)+(dphidphi(2)-dphidphi(1)-dQdphi(1)/(kD(2)/lam(2)*tanh(b(2)/lam(2))))/...
    B(1)*lam(1)/b(1)*sinh(b(1)/lam(1));
dhdphi(2)=dphidphi(2)+(dphidphi(1)-dphidphi(2)-dQdphi(1)/(kD(1)/lam(1)*tanh(b(1)/lam(1))))/...
    B(2)*lam(2)/b(2)*sinh(b(2)/lam(2));

% dq/dphi [(m/d)/m] is independent of recharge, just depends on fixed variables
dqdphi=b(1)/sum(b)/C(1)*(1-dhdphi(1))+b(2)/sum(b)/C(2)*(1-dhdphi(2)); % dq/dphi


%% Coordinats to compute heads

x1=0:0.1:b(1); % x in first compartment from left to right
x2=0:0.1:b(2); % x in second compartment from right to left

x =unique([x1 sum(b)-x2]);  % true real world x along cross section from center left

%% Computation of heads and mean heads and loop once to match q with
%  desired section-averaged seepage rate

NStep=2;  % we need exactly two steps to get the desired flux
hmean=NaN(1,2);
OutSteady=NaN(NStep,5);

%clc
NStep=5;

for i=1:NStep

    phi = (hLR(1,:)./csb+PHI./C)./(1./csb+1./C);% net fixed head in compartments

    QL=(N(1)*c(1)+phi(1)-hditch+...
    ((N(2)*c(2)+phi(2))-(N(1)*c(1)+phi(1)))*cosh(b(1)/lam(1))/B(1))/...
    (w+cosh(b(1)/lam(1))/(B(1)*kD(2)/lam(2)*tanh(b(2)/lam(2))));

    QR=(N(2)*c(2)+phi(2)-hditch+...
    ((N(1)*c(1)+phi(1))-(N(2)*c(2)+phi(2)))*cosh(b(2)/lam(2))/B(2))/...
    (w+cosh(b(2)/lam(2))/(B(2)*kD(1)/lam(1)*tanh(b(1)/lam(1))));
   
    Q=QL;
    
    % head as function of x1 and x2
     h1=N(1)*c(1)+((N(2)*c(2)+phi(2))-(N(1)*c(1)+phi(1))-Q/(kD(2)/lam(2)*tanh(b(2)/lam(2))))*cosh(x1/lam(1))/B(1)+phi(1);
     h2=N(2)*c(2)+((N(1)*c(1)+phi(1))-(N(2)*c(2)+phi(2))-Q/(kD(1)/lam(1)*tanh(b(1)/lam(1))))*cosh(x2/lam(2))/B(2)+phi(2);

 %    average head in both compartments
 %    hmean(1)=N(1)*c(1)+((N(2)*c(2)+phi(2))-(N(1)*c(1)+phi(1))-Q/(kD(2)/lam(2)*tanh(b(2)/lam(2))))*(lam(1)/b(1))*sinh(b(1)/lam(1))/B(1)+phi(1);
 %    hmean(2)=N(2)*c(2)+((N(1)*c(1)+phi(1))-(N(2)*c(2)+phi(2))-Q/(kD(1)/lam(1)*tanh(b(1)/lam(1))))*(lam(2)/b(2))*sinh(b(2)/lam(2))/B(2)+phi(2);
    
    hmean(1)=N(1)*c(1)+beta(1)*((N(2)*c(2)+phi(2))-(N(1)*c(1)+phi(1))-Q/(kD(2)/lam(2)*tanh(b(2)/lam(2))))+phi(1);
    hmean(2)=N(2)*c(2)+beta(2)*((N(1)*c(1)+phi(1))-(N(2)*c(2)+phi(2))-Q/(kD(1)/lam(1)*tanh(b(1)/lam(1))))+phi(2);

    % section-averaged seepage
    q=b(1)/sum(b)*(PHI-hmean(1))/C(1)+b(2)/sum(b)*(PHI-hmean(2))/C(2); % seepage
    

    hv=hmean(1)/(1+b(2)/b(1)*C(1)/C(2))+hmean(2)/(1+b(1)/b(2)*C(2)/C(1));
    cv=C(1)*C(2)*sum(b)/(b(1)*C(2)+b(2)*C(1));
 
    % show
    fprintf('PHI(%2d)=%10f  q=%10g hv=%10.2g cv=%10g (phi-hv)/cv=%10g  %12f QL=%12f QR=%12f hditch=%12f\n',...
         i,PHI,q,hv,cv,(PHI-hv)/cv,Q,QL,QR,hditch);

    % keep
    OutSteady(i,:)=[PHI q hv cv (PHI-hv)/cv];

    % update PHI for next loop
    PHI=PHI+(qsoll-q)/dqdphi;
    
end
PHI=OutSteady(end,1); % set PHI back to last used value

%% Plot results

figure; hold on; grid on; xlabel('x [m]'); ylabel('head [m]');
title(sprintf('wide ditch, q=%.3f m/d',q));
leg='';

plot(x1,h1,'r',sum(b)-x2,h2,'r'); leg{end+1}='h1'; leg{end+1}='h2';
plot([0 b(1) b(1) sum(b)],[hmean(1) hmean(1) hmean(2) hmean(2)],'g'); leg{end+1}='hmean';
plot([0 sum(b)],[PHI PHI],'k');                    leg{end+1}='PHI';

legend(leg)

x =[x1 fliplr(sum(b)-x2)];
hx=[h1 fliplr(h2)];

