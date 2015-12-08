function [GLG,GVG,GHG,HY]=getGXG(h,t,varargin)
%GETGXG computes GXG, i.e. average lowest, spring and highest groundwater
% level over hydrological years in the input.
% A hydrological year runs from April 1 through March 31
%
% USAGE:
%    [GLG,GVG,GHG,HY]=getGXG(h,t[,'plot',I])
%
% INPUT:
%  h = h(Nsec,Nt) where Nsec number of cross sections each of which
%                 representing a parcel between two parallel ditches
%                 Nt the number of times (days) for which we have h.
%  t = vector of times, (days) at th end of which we have this h
%  I = vector of section nrs or logical vector indicating which parcels to
%  plot.
%  options:
%  'plot' with possible arguments 'none' or false, 'time' or GXG
%         'time' plots GXG as function of time with GXG times indicaed in
%         the time series as markers.
%         'GXG'  plots GXGs next to each other
%
% used in ..
%   mflab/examples/mf2005/DutchTop/HSV
%   mflab/examples/mf2005/DutchTop/NPARK
%
% the heads are presumed ordered in the direction of the time vector.
% Therefore, a large number of time series can be handled at once
%
% TO 101023 151119

[Iplot,varargin] = getProp(varargin,'plot',[1 2]);

DV=datevec(t(:));  % [Y M D h min sec]

% hydrological years
startYr = DV(1,1);
if t(1)>datenum(startYr,4,1), startYr = startYr + 1; end

endYr   = DV(end,1) -1; % last hydrological year
if t(end)<datenum(endYr,3,31),  endYr = endYr-1; end

% run over hydrological years

hyYr = startYr:endYr; % hydrological startyr:endyr

for iy = numel(hyYr)-1:-1:1
    I = find(t>= datenum(hyYr,4,1) & t<= datenum(hyYr+1,3,31));
    [h_,I_] = sort(h(I),'ascend');
    HY(iy).Ighg  = I_(    1:  3);
    HY(iy).Iglg  = I_(end-2:end);
    HY(iy).glg   = h_(    1:  3);
    HY(iy).ghg   = h_(end-1:end);
    
    % Store the heads measured on 14th or 28th each month
    
    % There may be an issue if in a real 14 day measurement series the
    % 14th and 28th dates shift a bit due to weekend days. So to be robust
    % some play must be allowed for in selecting the date. In our case, we
    % have probably always daily data (because they are simulated on a
    % daily basis) and this issue does not play. But this must be realized
    % in general. TO 150419
    
    % Store all measurements on 14th or 28ty during current hydrologic year
    HY(i).J=find(t>=HY(i).t1 & t<HY(i).t2 & (DV(:,3)==14 | DV(:,3)==28));
    
    % Spring head:

    % Theo only takes april first
    % HY(i).K=find(DV(:,1)==years(i) & DV(:,2)==4 & DV(:,3)==1); % april first in hydrological years
   
    % However we now take the
    % the mean of march 14, march 28 and april 14
    % see http://www.grondwaterdynamiek.wur.nl/NL/Parameters/
    % Store respective indices in vector K
    HY(i).K = find(DV(:,1)==year & DV(:,2)==3 & DV(:,3)==14 | ...
                   DV(:,1)==year & DV(:,2)==3 & DV(:,3)==28 | ...
                   DV(:,1)==year & DV(:,2)==4 & DV(:,3)==14);
               
    % Sort the  heads on the 14th and 28th of current hydrological year
    % to allow taking the highest and lowest three of them
    % Also remember their locations [Is], to allow plotting them on their
    % correct time.
    [h_yr, Is] = sort(h(:,HY(i).J),2);   
    HY(i).GLG  = h_yr(:,1:3); % get the lowest three
    HY(i).GLGJ = Is(:,1:3);   % get their location in vector J
    
    HY(i).GHG  = h_yr(:,end-2:end); % get the highest three
    HY(i).GHGJ = Is(:,end-2:end);   % and their location in vector J
    
    HY(i).GVG=   h(:,HY(i).K);  % Get the spring head values
    HY(i).GVGJ=  HY(i).K;       % and the index in the overall series
end

GHG=mean([HY.GHG],2);
GLG=mean([HY.GLG],2);
GVG=mean([HY.GVG],2);


%% If no plot requested return
if isempty(varargin)
    return
end

%% Show data and points picked out for GXG computation

fprintf('plotting')
figure; hold on; xlabel(sprintf('%d - %d',years(1),years(end))); ylabel('head [m]'); grid on;
title('GLG(red), GVG(green) and GHG(blue)');
plot(t([1 end])',[GLG GLG],'r--');
plot(t([1 end])',[GVG GVG],'g--');
plot(t([1 end])',[GHG GHG],'b--');

for i=1:length(HY)
    plot(t(HY(i).J),  h(HY(i).J), 'bo');
    plot(t(HY(i).J(HY(i).GHGJ)),HY(i).GHG,'ro','markerFaceColor','r');
    plot(t(HY(i).K)            ,HY(i).GVG,'go','markerFaceColor','g');
    plot(t(HY(i).J(HY(i).GLGJ)),HY(i).GLG,'bo','markerFaceColor','b');
end
plot(t,h,'k');
datetick('x','yyyy');

legend('GLG','GVG','GHG','GLG','GVG','GHG','data');
