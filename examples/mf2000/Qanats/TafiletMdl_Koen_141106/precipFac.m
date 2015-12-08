%% This file computes the precipitation in the saemi arid areas of Morocco
% from 1100 as was derived by Till and Guiot (1990) based on year rings of
% cedar trees in different areas in Morocco.

% The graph was copied and digitized
% It is given in pixels as obtained from the overlaid graph drawn in PPT.
% Here it is reproduced and converted back to the original units.
% Then the averate precipitaion of the calibration period was added to turn
% the values into absolute precipitation values according to Till and
% Guiot.
% However, the precipitaiton in the calibration period for the semi-arid
% areas was 528 mm/a, whereas that in the Tafilalat over the twelve years
% that we have data of is 111 mm/d. Therefore to obtain a long-term series
% of precipitation for the Tafilalet, the values are downscaled to 111 mm
% over the calibration period.
% Finally a factor was computed as the actual value of yearly precipitation
% over that in the calibration period. This factor is used by
% multiplication with the measured precipitation over the 12 twelve year
% data period to obtain kind of historical data that reflect dry and wet
% multi year drought and wet periods in accordance with the data gathered
% by Till and Guiot.
% Given the fact that the correlation between the dry and wet periods was
% over 90% between the humid, semi-humid and semi-arid areas in Morocco, it
% is reasonable that a high correlation also applies for the arid area of
% the Tafilalet, and hence this method makes sense.

% TO 151011

TillGuiot = 528;  % mm/d average precip for semi-arid areas
Tafilet   = 111;  % mm/d Tafilalet for the 12 years that we have data for

load PrecMarocTillGuiot1990.txt

P = PrecMarocTillGuiot1990;

% figure;
% plot(P(:,1),P(:,2)); set(gca,'yDir','reverse');
t = 1100 + (P(:,1)-P(1,1)) / (P(end,1)-P(1,1)) * (1990-1100);

% rectify time, to make sure that time is monotnously increasing
dt     = diff(t);
I      = find(dt<=0);
t(I)   = (t(I)+t(I+1))/2;
t(I)   = t(I)-1;
t(I+1) = t(I)+2;

% Rescale from pixels digitized to factual precipitation data for the
% semi-arid areas in Morocco
% Notice that the values are relative to those in the calibraion period,
% that is the period 1924-1960, where their value of P is zero.
PTGmax    =  105; % max in mm/y in graph of Till and Guiot (1990) mm/y
PTGmin    = -175; % min in mm/y in graph of Till and Guiot (1990) mm/y
yPixMax =  201; % pixels
yPixMin =  309; % pixels

% Conversion of the pixels
PTillGuiot = PTGmin + (P(:,2)-yPixMin)/(yPixMax-yPixMin) * (PTGmax-PTGmin)  + TillGuiot; % [mm/y]

yStd = std(PTillGuiot);
PFac = PTillGuiot / TillGuiot; % precipiation factor

PTafilet = PFac * Tafilet;

figure; hold on;
plot(t,PTillGuiot); % precip according to Till & Guiot (1990) semi arid areas Morocco

plot(t([1 end]),TillGuiot -yStd * [0 0],'k');
plot(t([1 end]),TillGuiot -yStd * [1 1],'m');
plot(t([1 end]),TillGuiot -yStd * [2 2],'r');
plot(t([1 end]),TillGuiot +yStd * [1 1],'m');
plot(t([1 end]),TillGuiot +yStd * [2 2],'r');

title({'Pricip semi-arid regions Maroc 1100-1990, Till & Guiot (1990)'; ...
    'std relative to calbiration period (1924-1960)'});
xlabel('years'); ylabel('precip in mm/y');
legend('precip','+\sigma','+2\sigma','-\sigma','-2\sigma');

precipFac = [   t PTillGuiot PTafilet PFac; ...
             2000  TillGuiot  Tafilet 1.0];   % add for interpolation
t  = precipFac(:,1);
ty = (t(1):t(end,1))'; % time for all years regularly

precipFac = [ty, interp1(t,precipFac(:,2:end),ty)];

fprintf('%4.0f %4.0f %4.0f %8.3f\n',precipFac');

save precipFac precipFac