% Example Hemker and Maas (1987) figure 2 test case
% TO 131005

clear variables; close all;

basename='analVsMf'; % model basename

%% Input data case figure 2, Hemker and Maas (1987)

D=[1 1 1 1]'  ; kD  = [0 100 100 0]';     kh  = kD./D;
d=[1 1 1]';     c   = [100 100 100]';    kv  =  d./c;

Sf =[0   1   1   0]' * 1e-4;
St =[1.6 1.6 1.6  ]' * 1e-3;

% Switch these two lines on to attribute the storage of the aquitards
% to that of the layers to compare with ancient formulas that neglect
% storage inside aquitards
%Sf = Sf + [0; 0.5*(St(1:end-1)+St(2:end)); 0];
%St(:) = 1e-8;

sf=Sf./D;
st=St./d;


r=10; % piezometer distance used in graph by Maas/Hemker

%% model gird
xGr = [0 logspace(-1,4,50)];
yGr = [-0.5 0.5];

Nlay = numel(D)+numel(d);
dz   = zeros(Nlay,1); dz(1:2:2*numel(D)-1) = D; dz(2:2:2*numel(d))=d;
zGr = [0; -cumsum(dz)];

gr = gridObj(xGr,yGr,zGr,'AXIAL',true);

%% Model arrays
IBOUND = gr.const(1); IBOUND(:,:,[1 end])=-1; % only top and base fixed
STRTHD = gr.const(0);

Kh = zeros(size(dz)); Kh(1:2:2*numel(D)-1)=kh; Kh(2:2:2*numel(d))=kv;
Ss = zeros(size(dz)); Ss(1:2:2*numel(D)-1)=sf; Ss(2:2:2*numel(d))=st;

Kv = Kh;
Kh(2:2:2*numel(d))   =0;
Kv(1:2:2*numel(D)-1) =1e6;

HK    = gr.const(Kh);  
VK    = gr.const(Kv);
SS    = gr.const(Ss);

%% get wells to find extractions required for the analytic solution
well = wellObj(basename,'wells',gr,HK,'PER');

% Get the extractions for the analytic solution
for iw=numel(well):-1:1, Q(iw,1) = well(iw).Q(1); end

%% set up the fm2ct finite diference model (Matlab inplementation)

IH = zeros(gr.Nz,gr.Nx);
FH = zeros(gr.Nz,gr.Nx);
FQ = zeros(gr.Nz,gr.Nx); FQ([3 5],1)=Q;

[PERnams,PERvals] = getPeriods(basename);
T = PERvals(:,strmatchi('PERLEN',PERnams));
t = cumsum(T);

% run it to get the heads
[Phi ]= fdm2t(gr.xGr,gr.zGr(:),t,Kh,Kv,Ss,XS(IBOUND),IH,FQ,true);


%% Compute the drawdowns according to Theis and Hantush for one layer
sTheis1   = -Q(2)/(4*pi*kh(3)*D(3)) * theis(r^2*sf(3)./(4*kh(3)*D(3)*t));
sTheis12  = -sum(Q)/(4*pi*sum(kD(2:3))) * theis(r^2*(sum(sf(2:3))+st(2))./(4*sum(kD(2:3))*t));
lambda1   =  sqrt(kD(3)/(1/c(2)+1/c(3)));
sHantush1 = -Q(2)/(4*pi*kh(3)*D(3)) * hantush(r^2*sf(3)./(4*kh(3)*D(3)*t),r/lambda1);
lambda12  =  sqrt(sum(kD(2:3))/(1/c(1)+1/c(3)));
sHantush12 = -sum(Q)/(4*pi*sum(kD(2:3))) * hantush(r^2*(sum(sf(2:3))+st(2))./(4*sum(kD(2:3))*t),r/lambda12);


%% Here we generate as many layers within each aquitard as we like

layersPerAquitard = 6;

% collect the modle in modelObjects to allow subdividing layers
% automatically
Model(1) = modelObj(gr,'gridObj');
Model(2) = modelObj(HK,'3Dlay');
Model(3) = modelObj(VK,'3Dlay');
Model(4) = modelObj(SS,'3Dlay');
Model(5) = modelObj(STRTHD,'3Dlay');
Model(6) = modelObj(IBOUND,'3Dlay');
Model(7) = modelObj(HK,'3Dlay');

planes = 1:gr.Nlay+1; % we keep all planes (layer interfaces)

% define subdivisions per layer
subdivisions = [1 layersPerAquitard 1 layersPerAquitard 1 layersPerAquitard 1];

% subdivide the layers, do this for all elements in Model
Model = Model.changeLayers(planes,subdivisions);

% put the new model into the workspace, overwriting the old one
unpack

%% get the wells
well = wellObj(basename,'wells',gr,HK,'PER');

%% save what we need in mf_analyze
save underneath Q t St c Sf kD r Phi sTheis1 sTheis12 sHantush1 sHantush12 layersPerAquitard