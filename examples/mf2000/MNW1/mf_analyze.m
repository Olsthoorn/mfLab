%% mf_analyze example MNW1 Halford and Hansen (2002, p 15ff)
%  The documentation is unclear. This is as far as I could figure it out
%  and understood it. It has taken me a lot of time due to incomplete
%  documentation.
%  TO 110822

load('name'); load(basename); load underneath


H=readDat([basename,'.hds']);
B=readBud([basename,'.bgt']);
NT=length(B);

hrange=ContourRange(H,50);

%% Contour heads in layer 1 with MNW's, drains and chd (end situation)

figure('name','MNW1 example Halford and Hansen, 2002');

cmap = colormap;

iLAY = 1; SP=2;

for iLay = 1:2
    subplot(1,2,iLay,'nextplot','add','color',cmap(1,:));

    title(sprintf('Heads in layer %d, stress period %d',iLay,SP));
    xlabel('x [ft]'); ylabel('y [ft]');
    
    % Contour heads
    [c,hdl]=contourf(gr.xc,gr.yc,H(SP).values(:,:,iLay),hrange);
    clabel(c,hdl);
    colorbar;
    
    gr.plotGrid('edgecolor','w','edgealpha',0.15);

    % if WEL is on plot wells
    for iw=1:length(well)
        plot(well(iw).x,well(iw).y,'k.');
    end
    
    gr.plotBCN('DRN',DRN,iLay,SP,'b','linewidth',2);
    gr.plotBCN('CHD',CHD,iLay,SP,'r','linewidth',3);

    axis('equal'); axis('tight');

end

%% Print zone budgets for all zones

% Make a zone array using the well Nrs as zone numbers
ZONE=zeros(size(IBOUND));

ZONE(IBOUND==iCHD) = iCHD;
ZONE(IBOUND==iDRN) = iDRN;

for i=1:length(well), ZONE(well(i).idx)  = well(i).nr; end

%zonebudget(B,ZONE)

%% extract the MNW extractions and he heads over time direclty from the budget file

%Get the extraction of each of the wells
for iw = length(well):-1:1
    leg{iw} = well(iw).name;
    switch class(well)
        case 'wellObj', TYPE='WEL';
        case 'MNW1Obj', TYPE='MNW';
        case 'MNW2Obj', TYPE='MNW';
        otherwise
    end
    for it = length(B):-1:1
        QWel(iw,it) = sum(B(it).term{strmatchi(TYPE,B(it).label)}(ZONE==well(iw).nr));
        Head(iw,it) = mean(H(it).values(ZONE==well(iw).nr));
    end
end

time = [H.totim];

figure('name','MNW1 example Halford and Hansen, 2002');

subplot(2,1,1,'nextplot','add');
title('Multinode well limited extraction over time');
xlabel('time [d]'); ylabel('Discharge [ft3/d]');
plot(time,QWel);
legend(leg);

subplot(2,1,2,'nextplot','add');
title('Multinode well limited head  over time');
xlabel('time [d]'); ylabel('head [ft]');
plot(time,Head);
legend(leg)
