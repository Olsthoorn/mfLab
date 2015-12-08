%% Analyzing output of the model
% TO 120110
 
load('name.mat') % loads name.mat, which only contains the variable "basename"
load(basename);  % having retrieved the baasename value, we can load
load underneath;

H=readDat([basename,'','.hds']);  % use readDAT to read the heads  output file

%% Plot setup
figure; hold on;
title('red:Analytic, blue:Modflow, green:Fdm2c');
xlabel('x [m]');  ylabel('head [m]');
set(gca,'xscale','log');

%% Modflow heads in blue
for it=1:length(H)
    plot(gr.xm,H(it).values(:),'b');    
end

%% Compute analytic heads (Theis) and plot in red

t=[H.totim]';
u=S./(4*T*t)*gr.xm.^2;
s=well.Q/(4*pi*T)*expint(u);

plot(gr.xm,s,'r');

%% Compute the head using the matlab function fdm2ct.m and plot in green

FH=NaN(gr.Nlay,gr.Nx); FH(end)=0;
FQ=zeros(size(gr.xm));
FQ(1)=well.Q;

IH=zeros(size(FH));  % s(1,:);          % initial head

Phi = fdm2t(gr.xGr,gr.zGr(:),t,T,T,S,IBOUND,IH,FQ,'axial'); Phi=XS(Phi); % Theis numerically

plot(gr.xm,Phi,'g'); % Theis numerically  as function of time

%% Show that area are inccurate for low values of the drawdown

set(gca,'yscale','log');

%% Conclusion
fprintf('MODFLOW and fdm2ct give the same results but these are inaccurate for low\n');
fprintf('values of the drawdown. Number of interations is irrelevant\n');
