% Tutorial working with GoogleMap images
%
% SEE ALSO:
%   testURL testURL2 under mflab/examples/Visualization
%
% TO 110514

clear variables; close all;

% To work with GoogleMaps images in Matlab, you first have to generate an
% URL with the request. This URL requires the center of the picture, the
% zoom level and the size in pixels.
% For instance, to obtain a picture of a small island from Dubai World

center  = [ 25.253159, 55.178000];            % Dubai world island
zoom    = 16;           % zoom level must be between 0 and 21
pixsize = [600 400]';   % number of pixles must be between 1 and 640
format  = 'PNG';
maptype ='satellite';   % one of {'roadmap','satellite','terrain','hybrid'})

%%
% This is a way to define some paths that have been clicked and saved as kml files using GE.
% They will be plotted on the map.
kmlspec={ % kmlfilename, pathcolor, linewidth, fillcolor
    'A Dubai Wolrd Island.kml'                    ,'red'   ,1,'yellow';...
    'Reservoir Small Dubai World Island.kml'      ,'yellow',1,'yellow';...
    'Square600x600.kml'                           ,'yellow',1,'yellow'...
     };

%%
% We make each path a structure and connect it to the structure url, this way we can have
% as many paths as we like.

for j= size(kmlspec,1):-1:1
    path(j).color      =  kmlspec{j,2};
    path(j).weight     =  kmlspec{j,3};
    path(j).fillcolor  =  kmlspec{j,4};
    path(j).LATLON     =  kmlpath(kmlspec{j,1});
end

%% Get the googleMapObj
%
% Get the url that will trigger Google Maps to deliver the figure centered
% around the specified center, using the specified zoom level and having
% the size in zoomlevel-tiles coordinates (see howToManual under mflab/Doc/howtoManual)
% mf_GM2PIC yields the figure in the browser if no output argumeng its
% specified. If an output argument is specified, it yields just the URL
% that it would otherwise put in the webbrouwser's URL window.
% This is what we do here, because we just want the URL for now.

GM = googleMapObj(center,zoom,pixsize,maptype,format,'',path);

%%
%surf(GM.A);

%% Get size of figure in meters using tile pixel coordinates

[dx,dy]=GM.tileSize();

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

[xGr,yGr,zGr,xm,ym,ZM,DX,DY,DZ,NX,NY,NZ]=modelsize3(xGr,yGr,zGr); %% Coordinate housekeeping

[XM,YM]=meshgrid(xm,ym);  % full size matrices needed for HFB

Z=repmat(zGr,[NY,NX,1]);

axis equal
axis tight

plotgrid(xGr,yGr);

%% Plot the island and see if the coordinates are correct
%  We put these lines at the end to make sure they overlay the model grid
plot(island(   :,1),island(   :,2),'g','linewidth',1);
plot(reservoir(:,1),reservoir(:,2),'r','linewidth',1);

plot(square(   :,1),square(   :,2),'r','linewidth',1);

%% And so on
