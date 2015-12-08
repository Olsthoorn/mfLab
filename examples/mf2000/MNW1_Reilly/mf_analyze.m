load underneath
load name
load(basename);


%% UNDER CONSTRUCTION TO 110809

%% Analyzing model output

facealpha=0.3;

clr='w';

%% Reading unformattted files

H  =readDat([basename,'.hds']);
DDN=readDat([basename,'.ddn']);
C  =readMT3D('MT3D001.UCN');
%B=readBud([basename,'.bgt']);
%NT=length(B);
figure;
h=mf_3Dblock(gr.xGr,gr.yGr,gr.Z,C(3).values,20,10,5,facealpha,[0.8 0.8 0.8]);

view(3);

%% Figure

figure; hold on;
title('MNW1 example Reilly (See Konikow and Hornberger, (2006)'); xlabel('x [ft]'); ylabel('y [ft]');

hrange=ContourRange(H  (end).values(:,1:100,:),0.001);
drange=ContourRange(DDN(end).values(:,1:100,:),0.001);
crange=ContourRange(C,1);

[~,hdl]=contourf(gr.xm(gr.xm<250),XS(gr.zm),XS(H(end).values(end,gr.xm<250,:)),hrange);
set(get(hdl,'children'),'edgecolor','none');

%contourf(xm(xm<250),XS(zm),XS(DDN(end).values(end,xm<250,:)),drange);

for iw=1:numel(MNW)
    MNW(iw).plot2D('k',2)
end

%% Concentrations

% % THIS IS STILL UNDER CONSTRUCTION TO 10822
% figure; hold on
% %contourf(xm(xm<250),XS(zm),XS(C(end).values(end,xm<250,:)),crange);
% 
% contourf(xm,ym,C(end).values(:,:,1),crange)
% 
% contourf(xm(xm<25 & xm>-25),ym(ym<25),STCONC(ym<25,xm<25 & xm>-25,1),10); colorbar;

%%
figure; hold on;
contourf(gr.xm(gr.xm<50 & gr.xm>-50),gr.ym(gr.ym<50),C(4).values(gr.ym<50,gr.xm<50 & gr.xm>-50,2),crange); colorbar

