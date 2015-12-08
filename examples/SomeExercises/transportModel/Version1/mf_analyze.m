%% mf_analyze, used to analyze and visualize the output of the model

%% retrieve the basename of this model
load name      % retrieve basename stored in file name.mat
load(basename)  % get save model arrays

%% read the heads
H = readDat([basename '.HDS']);  % extension is .HDS see worksheet NAM

%% read the budget file
B = readBud([basename '.BGT']); % extension is .BGT see worksheet NAM

%% The heads are in H(it).values because there is only one time step
% H(it).values is the same as H.values

%% Ways to show the data

% showing a variable in a 3D grid using plotMesh
figure;  ax(1) = axes('nextplot','add');

% merging HK and VKCB to get a 3D array for both the layers and confining beds
% wihtout merging, you can just show the HK or the VKCB leaving non HK or
% non VKCB layers empty
var = NaN(gr.Ny,gr.Nx,gr.Nz);
var(:,:,gr.ITlay) = HK;
var(:,:,gr.ITcbd) = VKCB;

h = gr.plotMesh(ax(1),var,'title','conductivities','fontsize',14); view(3); % show it
%h = gr.plotMesh(ax(1),HK); view(3); % show it
%h = gr.plotMesh(ax(1),VKCB); view(3); % show it

colorbar;

%% plotting conductivities in a cross section along x-axis and along y-axis
gr.plotXSec(1,'figure','xsec','title','Conductivities along y-axis','fontsize',14,'all','lay',HK,'cbd',VKCB);
colorbar;

% plot the heads in the XSec figure
plot(gr.xm,XS(H(end).values(1,:,1)),'k','linewidth',2);
plot(gr.xm,XS(H(end).values(1,:,2)),'b','linewidth',2);

%% same thing along the y-axis
gr.plotYSec(1,'figure','ysec','title','Conductivities along x-axis','fontsize',14,'all','lay',HK,'cbd',VKCB);
colorbar;

% plot the heads in the XSec figure
plot(gr.ym,YS(H(end).values(:,1,1)),'k','linewidth',2);
plot(gr.ym,YS(H(end).values(:,1,2)),'b','linewidth',2);


%% Plot stream function
ix = 1;
B = mf_Psi(B,ix,'y',gr);  % compute stream function
prange = ContourRange(B,50,'PsiY',gr);

%plot stream function
contour(YS(gr.YP(:,1,:)),YS(gr.ZP(:,1,:)),B(1).PsiY,prange);


%% contouring is not useful in this case



%% Show the particles

%% Show particles in 3D
figure;   hold on;   view(3);   xlabel('x [m]');   ylabel('y [m]');   zlabel('z [m]');

gr.plotMesh('faceAlpha',0.15); % thin grey lines
pGrp.plot(); title('Particles starting points');

pGrp = pGrp.getEndPoints('mp6.endp');

%pGrp.dispEndPoints();

pGrp.endPointStatistics;

%pGrp.plotEndPoints(1,'ko');
pGrp.plotStartPoints(2,'b.');
pGrp.plotStartPoints(3,'r.');
pGrp = pGrp.getPathLines('mp6.plines');
pGrp.plotPath('b');
view(3);

