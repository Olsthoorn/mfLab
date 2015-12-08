%% Analyzing output of the model
% TO 120110
 
load('name.mat') % loads name.mat, which only contains the variable "basename"
load(basename);  % having retrieved the baasename value, we can load
load underneath;

H=readDat([basename,'.hds']);  % use readDAT to read the heads  output file
B=readBud([basename,'.bgt']);  % budget file
B=mf_Psi(B);

%% Plot setup
figure; hold on;
title('red:Analytic, blue:Modflow, green:Fdm2c');
xlabel('x [m]');  ylabel('head [m]');
set(gca,'xscale','log');

t=[H.totim];

%% Modflow heads in blue

figure; hold on; xlabel('x [m]'); ylabel('z [m]');

hrange=ContourRange(H,50);
prange=ContourRange(B,50,'Psi');

for it=1:length(H)
    if it==1
        [~,hdl1] = contour(gr.xc,gr.zc,XS(H(it).values),hrange,'r');
        [~,hdl2] = contour(gr.xp,gr.zp,B(it).Psi,prange,'b');
    else
        set(hdl1,'zdata',XS(H(it).values));
        set(hdl2,'zdata',B(it).Psi(2:end,:));
    end
end

%%
figure('name','Boulton drawdown'); hold on;
for it = 1:length(H)
    if it==1
        hdl = plot(gr.xm,XS(H(it).values)); set(gca,'xscale','log');
    else
        for j=1:numel(hdl)
            set(hdl(j),'ydata',XS(H(it).values(1,:,j)));
        end
    end
end

