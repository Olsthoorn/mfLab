%TESTHANTUSHNSCRIPT script that verifies the analytical mutilayer nantushn
%
% Example:
%    TestHantushnScript;
%
% Test analutical muli-layer solution by comparing it with outcomes of a
% numerical FD model that solves the same problem as given in the PhD
% thesis of Kick Hemker (2000), p Chapter 4, p 58 fig 3.
% This examaple is from the paper
%    Hemker and Maas (1987) in Journal of Hydrology
%
% The solution was first implemented by Fritz and Tijsen (2009) but
% I had to change the matrix function evaluation due to the change of
% funm in Matlab versions since 2007. It now uses eigen values and eigen
% vectors instead.
%
% See also hantushn stefest
%
% TO  090329

clear variables
close all

clr='bgrymck';
mrk='o+x*^v<>psh';

%% layers
layer={
'sdp1'    0  -10 1000 3e-3  
'wvp1'  -10  -30 2000 1e-3
'sdp2'  -30  -40 1500 5e-4
'wvp2'  -40  -60 1500 4e-4
'sdp3'  -60  -70 1000 3e-4
'wvp3'  -70  -90  500 1e-4
'sdp4'  -90 -100 4000 2e-4
'wvp4' -100 -120 2000 3e-4
'sdp5' -120 -130 20000 1e-3
};

%% Mesh and extraction
z=layer{1,2}:-5:layer{end,3};
r =logspace(1,log10(15000),50);
[r,z,rm,zm,dr,dz,Nr,Nz]=modelsize(r,z);

Q =[0, 10000, 0, 0];
t =[1e-3 1e-2 1e-1 1 10];

%% period
Period=[t(end) 50 1.25];

%% layer properties attribution

Kr=NaN(Nz,Nr);
Ss=NaN(Nz,Nr);

%% Aquitards
m=0;
c =zeros(ceil(size(layer,1)/2),1);
St=zeros(ceil(size(layer,1)/2),1);
for iLay=1:2:size(layer,1)
    m=m+1;
    c(m)=layer{iLay,4};
    St(m)=layer{iLay,5};
    D=layer{iLay,2}-layer{iLay,3};
    k=D/layer{iLay,4};
    ss=layer{iLay,5}/D;
    I=find(zm>layer{iLay,3} & zm<layer{iLay,2});
    Kr(I,:)=k;
    Ss(I,:)=ss;
end

%% Aquifers
m=0;
kD =zeros(floor(size(layer,1)/2),1);
Sf =zeros(floor(size(layer,1)/2),1);
for iLay=2:2:size(layer,1)
    m=m+1;
    kD(m)=layer{iLay,4};
    Sf(m)=layer{iLay,5};
    D=layer{iLay,2}-layer{iLay,3};
    k=layer{iLay,4}/D;
    ss=layer{iLay,5}/D;
    I=find(zm>layer{iLay,3} & zm<layer{iLay,2});
    Kr(I,:)=k;
    Ss(I,:)=ss;
end

%% Correction for boundary conditions
Kr(1,:)=Kr(1,:)./2;  % fixed head flow through half layer thickness only


%% FQ, 1 stress period
I=find(zm>layer{4,3} & zm<layer{4,2});
fq=Q(2)*dz(I)/sum(dz(I));
FQ=ones(length(I),5); FQ(:,3)=I; FQ(:,end)=fq;

%% FH, 1 stress period
FH=NaN(Nz,Nr); FH(1,:)=0;

%% IH initial head
IH=zeros(Nz,Nr);

%% Run Hantush
dd1=hantushn(Q,r,t,St,c,Sf,kD);

%% Run numerical model
[Phi,tnum,Qt]=fdm2tsp(r,z,Period,Kr,[],Ss,IH,FH,FQ,'radial');

dd2=interp3(rm,zm,tnum,Phi,rm,zm,t);  % interpolate to same times as analytical solution

%% Plot


IAquif=zeros(size(kD));  % first fined the numbers of the layers of the num model that correspond to the analytical aquifer numbers
m=0;
for i=2:2:floor(size(layer,1))
    m=m+1;
    IAquif(m)=find(zm<layer{i,2},1,'first');
end

for it=1:length(t)
    figure;
    semilogx(r,permute(dd1(:,it,:),[1 3 2]),'-'); hold on              % Plot head in analytical aquifers (lines)
    semilogx(rm,dd2(IAquif,:,it),'+');                                 % Plot head in numerical  aquifers ('+')
    title(sprintf('Hemker & Maas 90(1987)231-49 (analyt - vs fdm +) time = %g d',t(it)));
    xlabel('time [d]'); ylabel('drawdown [m]');
    set(gca,'ydir','reverse');
    grid on
end

% QED  we showed that the analytical and numerical model give the same
% result
% we alsd showed that the interpolation is efficient to match corresponding
% coordinates in space and time
% This implies that the analytical solution of Hantushn may now be used to
% verify models and for calibration
% TO 090329