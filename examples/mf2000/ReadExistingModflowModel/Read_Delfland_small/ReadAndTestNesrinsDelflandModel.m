% PH Candiate Nesrin worked on a groundwater model for the area around Delft, Netherlands,
% called the DelfLand model. She has a large and a small version. This is the smaller version.
% The model is given in MODFLOW input files. It is read-in using the mfLab function readMFLOW.
% After the model has been read, some checks are done, which were necessary to figure out her
% boundry conditions, and which resolved her model problem.
% This example is just to show how an existing model can be readin and then saved
% into a Matfile.
% Clearly, the matfile can be loaded by an mf_adapt file and thus be used to run
% the model through mfLab as usual.
%
% Notice the readMFLOW function
%
% TO 120323


%% First explore the heads
H=readDat('delft_14L_25m_tr.hed');

[Ny,Nx,Nz]=size(H(1).values);

t=[H.totim]; NPER=length(t);

%% The grid is not contained in the binary head file. It is generated here
%  externally using known values of dx, dy
dx=25; dy=25; xLL=0; yLL=0;

xGr=xLL+dx*(0 :  Nx);  xm=0.5*(xGr(1:end-1)+xGr(2:end));
yGr=yLL+dy*(Ny:-1:0)'; ym=0.5*(yGr(1:end-1)+yGr(2:end));

%% To explore the head easily: Rigorously make one 4D array (Ny,Nx,Nz,Nt) of head values
V=NaN(Ny,Nx,Nz,NPER);
for iPer=1:NPER
    V(:,:,:,iPer)=H(iPer).values;
end

%% Then plot the heads from different perspectives

iy=60; ix=70; iz=10; % choose a row, columen and layer to show head over time

%% time series for (iy,ix) all layers
figure;
plot(t,squeeze(V(iy,ix,:,:))); % time series of node ny,nx
title(sprintf('all layers at ix=%d iy=%d as function of time',ix,iy));
xlabel('time [d]'); ylabel('head [m]'); grid on

%% cross section one row, all layers at end of simulation
figure;
plot(xm,permute(V(iy,:,:,end),[3,2,1,4]))
title(sprintf('All layers at row %d at end of simulation t=%g',iy,t(end)))
xlabel('x [m]'); ylabel('head [m]'); grid on

%% map of heads for a single layer
figure;
[h,c]=contour(xm,ym,V(:,:,iz,end),-15:0.25:0); clabel(h,c);
title(sprintf('Contour of layer %d at end of simulation t=%g',iz,t(end)));
xlabel('x [m]'); ylabel('y [m]');grid on

%% Reading in the model into Modflow
% Notice that the model arrays is saved in model.mat at the end of readMLOW

readMFLOW

