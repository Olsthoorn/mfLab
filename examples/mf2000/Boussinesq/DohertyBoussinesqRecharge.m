% Example Boussinesq (flow on an sloped base)
% TO 100506

% This m files computes the steady state water table (head) in an unconfined aquifer
% with a partly steep slope under constant recharge. A Boussinesq problem.
% The main problem with this types of problems is converging difficulties
% due to cells alternatingly switching between wet and dry.
% Doherty (2001, Groundwater) has adapted MODFLOW to prevent dry cells
% altogether. He just gradually reduces the conductivity/transmissivity when
% the head is close to zero or below zero. This smooths the non-linearity
% that causes the alternating dry and wet cells.
% Because I don't have a copy of Doherty's MODFLOW (which seems to be
% included with the GMS user interface). However his method can easily
% be mimicked using the fdm3.m function, which implements a full fledge
% 3D steady-state finite difference model in MATLAB. It's confined so heads
% below the bottom of cells play no role. When put in a loop, the
% corrections necessary to simulate the unconfined flow can be done by
% adjusting the conductivities before the next loop cycle.
% The difference between MOFLOW and Matlab is, perhaps, that the MODFLOW
% PCG (and other) solvers solve iteratively while the Matlab solvers solves
% in its entirety without referring to iterations and hence to inner and
% outer loops. This is why the MODFLOW outerloop adjustments are done after
% each run of fdm3.m within an explicit loop.
%

clear variables;

%% The model name. Every model has a "basename" and all associated files

clr='brgkmc';
k =2;   % conductivity (uniform)
dx=1;   % grid cell width (uniform)
z0=0;   % z at xm(1)
D=10;   % dummy aquifer thickness (because unconfined)
slope =-1/2; % dB/dz, inclination with B the bottom of the aquifer
Prec  = [0.001 0.01 0.1];  % precipitation
%kstar=k;  % k*cos(theta)^2;

%% Specify grid line coordinates
xGr=0:dx:1000; % [m]
yGr=[-0.5 0.5];  % [m]
[xGr,yGr,xm,ym,DX,DY,Nx,Ny]=modelsize(xGr,yGr);

%% Aquifer bottom, top and center elevation

A=40; % set A=0 for constant slope
zB=z0+slope*(xm-xm(1))+A*sin(2*pi*(xm-xm(1))/(xm(end)-xm(1))*2);  % Aquifer bottom
zT=ones(size(xm))*D;  Nz=1;  % in face arbitrary top because unconfined

Z=NaN(Ny,Nx,Nz+1);
Z(:,:,1)    =zT;
Z(:,:,end)  =zB;

Dz=-diff(Z,1,3);

ZB=Z(:,:,2:end);                    % make sure the second aquifer is always filled but has small k

zm=0.5*(Z(:,:,1:end-1)+Z(:,:,2:end));

%% Arrays
IBOUND=ones(Ny,Nx,Nz);  IBOUND(:,[ end],1)=-1; % only the first layer

HK    =ones(Ny,Nx,Nz)*k; % hor conductivity      [m/d]

%% We set the starting heads equal to the normal depth
STRTHD=Z(:,:,1:end-1);                  % in fact irrelevant because steady state
STRTHD(IBOUND==-1)=ZB(IBOUND==-1)+D;    % right hand fixed head boundary

%% Try the same with fdm3
leg=cell(2*length(Prec),1);

figure; ax1=axes; hold on; plot(xm,squeeze(Z(:,:,end)),'color','k','linewidth',1);
figure; ax2=axes; hold on; plot(xm,squeeze(Z(:,:,end)),'color','k','linewidth',1);

m=0;  % counter of legend test
m=m+1; leg{m}='bottom of aquifer';

for iP=1:length(Prec)
    figure;
    FH=STRTHD; FH(IBOUND~=-1)=NaN;
    FQ=zeros(size(FH));
    FQ(:,:,1)=DX*Prec(iP);
    
    %% parameters to mitigate non-linearity of dry cells acc. Doherty 2001    
    alfa=5; % to fix d2 in terms of d0
    d0=2.5; % see this as capillary zone thickness
    d2=alfa*d0;           % d2 is chose in terms of d0
    d1=alfa*d0/(alfa-1);  % this fixes d1

    plot(xm,squeeze(Z(:,:,end)),'color','k','linewidth',1); hold on % plot base

    % Dohorty (2001 GW) rewetting method
    %% initialize
    N=1000;
    T=zeros(size(HK));
    B=min(Dz,(STRTHD-Z(:,:,2:end)));
    T(B> 0)=HK(B> 0)*d0.*exp(-B(B> 0)/d1)+HK(B>0).*B(B>0);
    T(B<=0)=HK(B<=0)*d0.*exp( B(B<=0)/d2);
    K=T./Dz;
    PhiOld=STRTHD;
    %% Run
    for iter=1:N;
        [Phi,Q,Qx]=fdm3(xGr,yGr,Z,K,K,K,FH,FQ);
        d2=alfa*d0;
        d1=alfa*d0/(alfa-1);
        Kold=K;
        B=min(Dz,(Phi-Z(:,:,2:end)));
        T(B> 0)=HK(B> 0)*d0.*exp(-B(B> 0)/d1)+HK(B>0).*B(B>0);
        T(B<=0)=HK(B<=0)*d0.*exp( B(B<=0)/d2);
        K=0.5*Kold+0.5*T./Dz;
        dPhi=max(abs(PhiOld(:)-Phi(:)));          fprintf('iteration %d, dPhiMax=%g m/d',iter,dPhi);
        I=find(dPhi==max(abs(PhiOld(:)-Phi(:)))); fprintf('Position of max Phi change %d, length(I)=%d\n',I(1),length(I));
        if dPhi<0.001, break; end % ready, if so jump out

        % renew figure ever 50 iterations to prevent mess/clutter
        if rem(iter,50)==0
            hold off
            plot(xm,Z(:,:,end),'color','k','linewidth',1); hold on
        end
        
        plot(xm,squeeze(Phi(:,:,  1)),'b'); hold on
        drawnow;
        PhiOld=Phi;
    end
    close(gcf);
    %% when done, switch to ax1 and plot the new final line to it
    axes(ax1); plot(xm,squeeze(Phi),clr(iP)),'--'; hold on
    axes(ax2); plot(xm,squeeze(Z(:,:,end)+T./HK),'color',clr(iP),'linewidth',1);
    m=m+1; leg{m}=sprintf('Prec=%g m/d',Prec(iP));
end
leg=leg(~cellfun('isempty',leg));

%% also for axis ax1, place labels and title
axes(ax1);
legend(leg);
xlabel('x [m]'); ylabel('elevation, head [m]'); grid on;
title(sprintf('head along slope with recharge [%g %g %g] m/d',Prec));

axes(ax2);
legend(leg);
xlabel('x [m]'); ylabel('elevation, head [m]'); grid on;
title(sprintf('head along slope with recharge [%g %g %g] m/d',Prec));
