% MF_ADAPT: Example use of HFB Package and set_HFB function
%
% USAGE:
%   mf_adapt  % to juse see of it works
%   mf_setop  % in invoke it and run the models
%
% SEE ALSO:
%   testURL testURL2 under mflab/examples/Visualization
%
% TO 110514

clear variables; close all;

basename='FSSE-HFB';  % example that uses this piece of code

%% Coordinates of an island from Dubai World
url.center = [ 25.253159, 55.178000];            % Dubai world island
url.zoom   = 16;           % between 0 and 21
url.size   = [600 400]';   % each between 1 and 640

url.maptype='satellite'; % one of

%% Some paths clicked and saved as kml files from Goolge Earth, 
kmlspec={
    'A Dubai Wolrd Island.kml'                    ,'red'   ,1,'yellow';...
    'Reservoir Small Dubai World Island.kml'      ,'yellow',1,'yellow';...
    'Square600x600.kml'                           ,'yellow',1,'yellow'...
     };

%% We make each path a struction to allow adding some extra info like color
for j=1:size(kmlspec,1)
    url.path(j).color      =  kmlspec{j,2};
    url.path(j).weight     =  kmlspec{j,3};
 %   url.path(j).fillcolor  =  kmlspec{j,4};
 
 %% kmpath extracts the path coordinates from the kml path file
    url.path(j).LATLON     =  kmlpath(kmlspec{j,1});
end

% Get the url that will trigger Google Maps to deliver the figure centered
% around the specified center, using the specified zoom level and having
% the size in zoomlevel-tiles coordinates (see howToManual under mflab/Doc/howtoManual)
% mf_GM2PIC yields the figure in the browser if no output argumeng its
% specified. If an output argument is specified, it yields just the URL
% that it would otherwise put in the webbrouwser's URL window.
% This is what we do here, because we just want the URL for now.

URL=mf_GM2PIC(url);                % URL to get the GM picture

%% Get size of figure in meters using tile pixel coordinates

[dx,dy]=mf_GM_tile_size(url.zoom,url.center(1));

xLim=dx*(url.size(1)/256)*[-0.5 0.5];  % every GM tile is 256x256 pixels
yLim=dy*(url.size(2)/256)*[-0.5 0.5];  % of known absolute size

%% Here we use the URL twice. First to get the figure into Matlab
% Then to get its information, becuase we need the colormap later on
A    = imread(URL);                 % get it directly into matlab
info = imfinfo(URL);               % get associated image info

%% Here the fiture is plotted, using the coordianates in m centered around
%% the center of the figure.
image(xLim,yLim([2 1]),A); colormap(info.Colormap); % Plot it using the colormap from the info

hold on;
set(gca,'ydir','normal');
xlabel('x [m relative to tile center]');
ylabel('y [m relative to tile center]');
title('A small Dubai World Island, suitable for subsurface water revervoirs?');

%% Now lets use the GE path by converting its Lat Lon coordinates into x,y
%% coordinates relative to the point url.center.
island=   mf_GMLL2XY(url.path(1).LATLON,url.center);
reservoir=mf_GMLL2XY(url.path(2).LATLON,url.center);

%% We may add a square of 600 by 600 m to allow verification with Google
%% Earth's ruler

square   =mf_GMLL2XY(url.path(3).LATLON,url.center);

%% Here we show how this can be used in your model

%% Grid

xGr=[-400:20:500 -200:10:400]; 
yGr=[-350:20:350 -100:10:200];
zGr=[3 -5:-2.5:-50];  % Plane elevation vector

gr = gridObj(xGr,yGr,zGr); %% Coordinate housekeeping

axis equal
axis tight

gr.plotGrid('c');

%% Plot the island and see if the coordinates are correct
%  We put these lines at the end to make sure they overlay the model grid
plot(island(   :,1),island(   :,2),'g','linewidth',1);
plot(reservoir(:,1),reservoir(:,2),'r','linewidth',1);

plot(square(   :,1),square(   :,2),'r','linewidth',1);

%% And so on
