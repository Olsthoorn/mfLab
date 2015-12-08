%% Analyzing output of the model
% TO 121126

clear; close all

load('name'); load(basename); load('underneath');

%% ====== REQUEST ANIMATE OBJ ======
animate = animateObj(basename,{'salinity','budget'});

if exist('MNW','var')
    well = MNW.setCout(animate.UCN{:});  % add computed conc for all species simultaneously
else
    well = well.setCout(animate.UCN{:});  % add computed conc for all species simultaneously
end

%% Animate concentration
figpos = get(0,'screensize'); figpos = [figpos([1 2]) 2/3*figpos([3 4])];

%% show temperature curves and conentration of all wells
well = animate.cOut(well);

%% animate conc in the plane
iLay = 1; iComp=1;
animate.concXY(gr,iComp,iLay,well,[],'mesh','on');

%% animate conc in the plane
iRow = well(1).iy(1); iComp=1;
animate.concXS(gr,iRow,well)
%% Check water balance for all wells
[WB1,time] = well.massBudget(1,[Inf Inf]);

%% Total extraction
clear w;

if strcmp(class(well),'wellObj'), lbl = 'WELLS'; else lbl='MNW'; end

B = readBud([basename '.BGT']);
for it=length(B):-1:1
    try
        imnw = strmatchi(lbl,B(it).label);
    catch ME
    end
    for iw=numel(well):-1:1
        w(iw).Q(it,:) = B(it).term{imnw}(well(iw).idx);
    end
end

for it=length(animate.UCN{iComp}):-1:1
    for iw=numel(well):-1:1
        w(iw).C(it,:) = animate.UCN{iComp}(it).values(well(iw).idx);
    end
end
