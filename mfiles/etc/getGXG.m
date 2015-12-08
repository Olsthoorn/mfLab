function [GLG,GVG,GHG,HY]=getGXG(h,t,varargin)
%GETGXG -- compute GLG GVG and GHG for the periodgiven
%
% USAGE:
%         [GLG,GVG,GHG]=getGXG(h,t[,'plot',I])
%
% h = heads(Nsection,Nt)
% t = time vector (1,Nt)
% 'plot', I  label value pair indicating which sections to plot vs time
%
% Compute the "average lowest head" (GLG), "spring head" (GVG) and "average
% "highest head" (GHG) from the given heads over the given time span.
% The heads are presumed to be ordered in the direction of the time vector.
% Therefore, a large number of time series can be handled at once
% plotMode=0 or omitted altogether, don't plot
% plotMode=1 plot heads over time with points that were used to compute CxG
% plotMode=2 plot GxG for all sections with section number on x-axis
%
% TO 101023 JB 101119 TO 101119

[Iplot, ~ ] = getProp(varargin,'plot',[]);


DV      =datevec(t(:));  % [Y M D h min sec]
startYr = DV(1,1);
endYr   = DV(end,1);

apr01startYr = datenum(startYr, 4, 1);

if t(1)>=apr01startYr, startYr = startYr + 1; end

endYr = endYr - 1;
years = startYr:endYr; % these are also the hydrological years

%% Hydrological years 1-4 tot 31-3

for i = numel(years):-1:1
    year = years(i);
    HY(i).t1 = datenum(year  ,4, 1); % start hydrological year
    HY(i).t2 = datenum(year+1,3,31); % end   hydrological year
    
    %% Get indices pointing to values inside the current hydrological year
    % ItYr are the indices in the original time series array of values that
    % fall within the curent hydrological year
    HY(i).ItYr  = find( t>=HY(i).t1 & t<=HY(i).t2 )';
    
    %% Store the heads measured on 14th or 28th each month
    
    % There may be an issue if in a real 14 day measurement series the
    % 14th and 28th dates shift a bit due to weekend days. So to be robust
    % some play must be allowed in selecting the date. In our case, we
    % have probably always daily data and this issue does not play. But
    % this must be realized in general. TO 150419
    
    % Store indices ponting at 14th or 28ty during current hydrologic year
    HY(i).It2wk = find( t>=HY(i).t1 & t<HY(i).t2 & (DV(:,3)==14 | DV(:,3)==28) )';
    
    %% Spring head:
    % Here defined as the mean of march 14, march 28 and april 14
    % in the same year (not hydrological year)
    % see http://www.grondwaterdynamiek.wur.nl/NL/Parameters/
    % Store respective indices in vector K
    HY(i).ItSpr = find(DV(:,1)==year & DV(:,2)==3 & DV(:,3)==14 | ...
                       DV(:,1)==year & DV(:,2)==3 & DV(:,3)==28 | ...
                       DV(:,1)==year & DV(:,2)==4 & DV(:,3)==14)';
               
    % Sort the  heads on the 14th and 28th of current hydrological year
    % to allow taking the highest and lowest three of them
    % Also remember their locations [Is], to allow plotting them on their
    % correct time.
    [h2wkSorted, Ist] = sort(h(:,HY(i).It2wk),2);  
    
    HY(i).GLG   = h2wkSorted(    :,1:3); % get the lowest (i.e. first) three values
    HY(i).GLGit = HY(i).It2wk(Ist(:,1:3));  % get their location in vector It2wk
    
    HY(i).GHG   = h2wkSorted( :,end-2:end); % get the highest three
    HY(i).GHGit = HY(i).It2wk(Ist(:,end-2:end));   % and their location in vector It2wk
    
    HY(i).GVG   = h(:,HY(i).ItSpr);  % Get the spring head values
    HY(i).GVGit=  ones(size(h,1),1)*HY(i).ItSpr;  % and the index in the overall series
end

% Compute the GXG over all sections from those of the hydrological years
GHG=mean([HY.GHG],2);  % GHG(Nsections,1)
GLG=mean([HY.GLG],2);
GVG=mean([HY.GVG],2);

%% Show data and points picked out for GXG computation  as function of time
if  isempty(Iplot), return; end
    
defaults = {'nextPlot','add','xGrid','on','yGrid','on'};

fsz = 12;

figure('position',screenPos(0.9));

ax = axes(defaults{:},'fontsize',fsz);

xlabel(ax,sprintf('%d - %d',years(1),years(end)+1));
ylabel(ax,'head [m]');
title( ax,...
    {sprintf('Heads for %d time series with GLG(red), GVG(green) and GHG(blue)',numel(Iplot));...
    'grey=daily, black o = 2-weekly, magenta = HYyear boundary'},'fontSize',fsz);

%% plot GHG, GLG and GVG als horizontal lines
plot(ax, t([1 end])',[GHG(Iplot) GHG(Iplot)]','b--');
plot(ax, t([1 end])',[GVG(Iplot) GVG(Iplot)]','g--');
plot(ax, t([1 end])',[GLG(Iplot) GLG(Iplot)]','r--');

plot(t,h(Iplot,:),'color',grey,'lineWidth',0.1);

%% plot the data points on which the GHG, GLG and GVG are based
for i=1:numel(HY)        
    plot(ax, t(HY(i).It2wk),h(Iplot,HY(i).It2wk), 'ko','markerSize',9); 

    for ip = 1:numel(Iplot)
        plot(ax, t(HY(i).GHGit(Iplot(ip),:)),HY(i).GHG(Iplot(ip),:),'bo','markerFaceColor','b');
        plot(ax, t(HY(i).GVGit(Iplot(ip),:)),HY(i).GVG(Iplot(ip),:),'go','markerFaceColor','g');
        plot(ax, t(HY(i).GLGit(Iplot(ip),:)),HY(i).GLG(Iplot(ip),:),'ro','markerFaceColor','r');        
    end
end

datetick();

% plot boundaries of hydrological years
yLim = get(ax,'ylim');
for i=1:numel(HY)
    plot(HY(i).t1([1 1]),yLim,'m');
    plot(HY(i).t2([1 1]),yLim,'m');
end

