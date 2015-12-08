% Example Hemker and Maas (1987) figure 3 test case
% TO 090806 091129

clear variables; close all;

basename='analVsMf'; % model basename

D=[1 1 1 1 1 1]'  ; T  = [0 2000 1500  500 2000 0]';     kh  = T./D;
d=[1 1 1 1 1]';     c  = [1000 1500 1000 4000 20000]';    kv  =  d./c;

Sf =[ 0  1  4  1  3  0]' * 1e-3;  sf=Sf./D;
St =[ 3  5  3  2  1 ]'   * 1e-3;  st=St./d;

xGr = [0 logspace(-1,5,50)];
yGr = [-0.5 0.5];

Nlay = numel(D)+numel(d);
dz   = zeros(Nlay,1); dz(1:2:2*numel(D)-1) = D; dz(2:2:2*numel(d))=d;
zGr = [0; -cumsum(dz)];

gr = gridObj(xGr,yGr,zGr,'AXIAL',true);

IBOUND = gr.const(1); IBOUND(:,:,[1 end])=-1;
STRTHD = gr.const(0);

Kh = zeros(size(dz)); Kh(1:2:2*numel(D)-1)=kh; Kh(2:2:2*numel(d))=kv;
Ss = zeros(size(dz)); Ss(1:2:2*numel(D)-1)=sf; Ss(2:2:2*numel(d))=st;

Kv = Kh;
Kh(2:2:2*numel(d))   =0;
Kv(1:2:2*numel(D)-1) =1e6;

HK    = gr.const(Kh);  
VK    = gr.const(Kv);
SS    = gr.const(Ss);

%% Here we generate as many layers within each aquitard as we like

layersPerAquitard = 10;

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
subdivisions = [1 layersPerAquitard ...
                1 layersPerAquitard ...
                1 layersPerAquitard ...
                1 layersPerAquitard ...
                1 layersPerAquitard ...
                1 ];

% subdivide the layers, do this for all elements in Model
Model = Model.changeLayers(planes,subdivisions);

% put the new model into the workspace, overwriting the old one
unpack


well = wellObj(basename,'wells',gr,HK,'PER');

for iw=numel(well):-1:1, Q(iw,1) = well(iw).Q(1); end


[PERnams,PERvals] = getPeriods(basename);
Dt = PERvals(:,strmatchi('PERLEN',PERnams));
t = cumsum(Dt);

save underneath Q t St c Sf T layersPerAquitard