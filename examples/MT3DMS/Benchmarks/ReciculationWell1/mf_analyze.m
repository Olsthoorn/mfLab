%% Analyzing output of the model
% TO 121126

clear; close all

load('name'); load(basename); load('underneath');

%% ====== REQUEST ANIMATE OBJ ======
animate = animateObj(basename,{'salinity','budget'});

if ~iswell, well=MNW; end

well = well.setCout(animate.UCN{:});  % add computed conc for all species simultaneously

%% Animate concentration
figpos = get(0,'screensize'); figpos = [figpos([1 2]) 2/3*figpos([3 4])];

%% show temperature curves and conentration of all wells
well = animate.cOut(well);

%% animate conc in the plane
iLay = 1; iComp=1;
animate.concXY(gr,iComp,iLay,well,[],'mesh','on');

%% Check water balance for all wells
[WB1,time] = well.massBudget(1,[Inf Inf]);

