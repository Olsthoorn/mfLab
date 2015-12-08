
%% generate testset
% 10 101113 110102

%% After the original data have been read in, change the data into a test

di=50;
if 1
    load tne
    tne(:,2)=tne(:,2)-tne(:,3); tne(:,3)=0;
    n=  1;    tne((n-1)*di+1:n*di,2)=  0.01;
    n=n+1;    tne((n-1)*di+1:n*di,2)=  0.0;
    n=n+1;    tne((n-1)*di+1:n*di,2)= -0.01;
    n=n+1;    tne((n-1)*di+1:n*di,2)=  0;

    tne(:,1)=(tne(:,1)-tne(1,1))/5;
    tne(:,2)=0.01;
end

if 1
   tne=tne(1:200,:);  % instead of 3288
   
%    for iP=1:length(P)
    clear P;
    C_DEK_MIN=1;

    N        =   0.01;  % recharge
    P.q      =    0;    % -0.003; % upward seepage

    P.h_mean   = 0; -0.6;
    P.h_winter = 0; -0.6;
    P.h_summer = 0; -0.6;
    P.phi    = -1.0; % head in regional aquifer (if used)


    P.b   = 50 ;  % width of ditch and width of parcel outside ditch
    P.D1  = 10;   %2;
    P.D2  = 10;   %50;

    P.AHN = 0.4;
    P.zdr = P.AHN;     % [NAP] default drain elevation    
    P.z0  = P.h_mean;
    P.z1  = P.h_mean-P.D1;  % [NAP] default bottom of top aquifer
    P.z2  = P.z1    -P.D2; % [NAP] default bottom of second, regional aquifer'

    mu    = 0.1;
    ss    = 1e-5;
  
    P.sy1 = mu;      % [ - ] specific yield of first layer
    P.S1  = mu;      % [ - ] confined storage coefficient (MF2005 STORAGECOEFFICIENT option)
    P.ss1 = ss;      % [1/m] specific storage coefficient
    P.sy2 = mu;      % [ - ] specific yield of first layer
    P.S2  = mu; % ss*P.D2; % [ - ] elastic storage coefficient layer 2 (MF2005 STORAGECOEFFICIENT option)
    P.ss2 = ss;      % [1/m] default elastic storativity, all layers

    P.hk1 = 10;
    P.hk2 = 10;  % 30;
    P.c   = 100;
    P.vk1 = 0.5 * P.D1/P.c; % [ - ] vert anisotropy first layer, requires layvka to be 1 in the LAY worksheet
    P.vk2 = P.hk2;           % [ - ] default vertical cond regional aquifer, require layvka=1 in the LAY worksheet

    P.cdr = 0.1;          % [ d ] default drain resistance

    %% Ditch properties
    P.dw   =   2;       % default ditch width
    P.cdb  =   1;    % 1;       % ditch default entry resistance value
    P.dd   = 2.5;       % [ m ]default ditch depth
    P.zdbot= P.h_mean-P.dd;
    if P.zdbot>P.z1
        P.Omega1 = min(P.dw/2+P.dd,P.h_mean-P.z1);
        P.Omega2 = 0;
        P.w1     = P.cdb*P.D1/P.Omega1;  %+2/(pi*sqrt(P.hk1*P.vk1))*log(P.D1/P.Omega1*sqrt(P.hk1/P.vk1));
        P.w2     = 1e6;  %+2/(pi*sqrt(P.hk2*P.vk2))*log(P.D2/P.Omega2*sqrt(P.hk2/P.vk2));
  
        P.vk_ditch = min(P.vk2, 0.5*(P.zdbot-P.z1)/P.D1 / P.c); %+ 2/(pi*sqrt(P.hk1*P.hk2))*log(2*P.D2/P.dw *sqrt(P.hk1/P.vk2));
        
    else
        P.Omega1= P.h_mean - P.z1;
        P.Omega2= P.dw/2+(P.z1-P.zdbot);
        P.w1    = P.cdb;                % because Omega1 = total depth of first aquifer
        P.w2    = P.cdb*P.D2/P.Omega2 +2/(pi*sqrt(P.hk2*P.vk2))*log(P.D2/P.Omega2*sqrt(P.hk2/P.vk2));
        
        P.vk_ditch = P.vk2;
    end

    DV=datevec(tne(:,1)); MONTH=DV(:,2);
    P.hLR=P.h_winter*ones(size(tne(:,1)));
    P.hLR(MONTH>=4 & MONTH<=9)=P.h_summer;
    
    %% Test slootpeil
    P.hLR=ones(size(P.hLR))*P.h_mean; % P.hLR(75:end)=P.hLR(75:end)+0.25; P.hLR(125:end)=P.hLR(125:end)-0.5;
    
    P.hLR=0*P.hLR;  % set hLR op nul

%    end
end