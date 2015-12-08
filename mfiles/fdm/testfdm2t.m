%TESTFDM2T test of numerical transient fdm model fdm2t by comparing with
%  drawdown according to the analytical Theis drawdown
%  Radial symmetric flow
%
% EXAMPLE
%    testfdm2t
%
% ToDo: unify these svn
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
%  TO 101211

%% Setup of problem
rW=0.1;      % well radius
k=10;        % conductivity
Ss=0.0001;   % specific storativity
Sy=0.1;      % specific yield
Qw=-2400;    % extracion

t=logspace(-3,2,21);       % times
r=logspace(log10(rW),4,41);% distances
z=[0, -20];                % top and bottom of aquifer

[r,z,rm,zm,dr,dz,Nr,Nz]=modelsize(r,z);	% grid coordinate housekeeping

S      = Ss*ones(Nz,Nr);   % specific storage in all cells (if more than 1 layer)
S(1,:) = Sy/dz(1);	       % overwrite top of model with specific yield
K      = ones( Nz,Nr).*k;  % conductivity field
FH     = NaN(  Nz,Nr);     % fixed head (locations), no fixed heads
Q      =-2400;             % well extraction
FQ     = zeros(Nz,Nr); FQ(:,1)=Q; % extraction field

kD=k*sum(dz);       % transmissivity field

IH=Q/(4*pi*kD)*expint(rm.^2*Sy/(4*kD*t(1))); % initial head (analytic)

%% Run transient finite difference model
[Phi,Qt,Qr,Qz,Qs]=fdm2t(r,z,t,K,K,S,IH,FH,FQ,'radial'); % fdm model

%% plot and compute analytical drawdown

figure; hold on;	% start visualisation for it=2:length(t)
for it=1:length(t)
    semilogx(rm,Phi(1,:,it),'ro');                                % numerical
    semilogx(rm,Q/(4*pi*kD)*expint(rm.^2*Sy/(4*kD*t(it))),'-');  % Theis
end

%% round off
set(gca,'xscale','log','xgrid','on','ygrid','on');
xlabel('time [d]'); ylabel('drawdown [m]');
title('numerical and analytical computation of Theis drawdown');
