% TEST_GM2PNG: EXAMPLE of how to use mf_GM2PNG to retrieve figure from Google Maps
% and put it on a figure with X Y coordinates relative to center
%
% USAGE:
%    test_GM_PIC_LL
%
% TO 110501 110514

clear variables; close all

%% We start making a struct holding the parameter values that Goolge Maps
%  requires as building blocks of a request for a map.
%  we need center, zoomlevel, size (in pix) and maptype

%    url.center='Brooklyn+Bridge,New+York,NY';
%    url.center='Dubai,AE';    % as name
%    url.center='The+Palm+Jebel+Ali,Dubai,AE';    % as name
url.center=[25.003395,54.986850];            % alternative

url.zoom=12;           % between 0 and 21
url.size=[512 512]';   % each between 1 and 640
url.maptype='roadmap'; % one of

%% Add as many markers to the struct as you like
i=0;

i=i+1;
url.markers(i).color='blue';
url.markers(i).label='S';
url.markers(i).LATLON=[25.002696,54.988728];

i=i+1;
url.markers(i).color='green';
url.markers(i).label='G';
url.markers(i).LATLON=[ 25.030855, 54.971621];

i=i+1;
url.markers(i).color='red';
url.markers(i).label='C';
url.markers(i).LATLON=[ 24.982417,54.970744];

%% Add as many paths to your struct as you like
url.path.color='blue';         % path color
url.path.weight=5;             % path line thickness
url.path.fillcolor='yellow';   % color to fill interior of path polygon

% actual path, as Lat Lon in decimal degrees
% you can use kmlpath('kmlpathfile.kml')  instead, where 'kmlpathfile.kml
% is the path clicked in Google Earth and saved as kml file.
url.path.LATLON=[
    24.968539  55.012772
    24.971400  55.009300
    24.973921  55.007122
    24.977864  55.005253
    24.988782  54.998962
    25.002550  54.990773
    25.014616  54.983683
    25.017455  54.981950
    25.030233  54.975013
    25.027804  54.969686
    25.022242  54.961701
    25.010288  54.955775
    24.999767  54.954395
    24.994937  54.955820
    24.988834  54.959602
    24.982031  54.967933
    ];
 
%    URL=mf_GM2PNG(url); % with output argument to just get the URL
[xW yS xE yN]=mf_GM_PIC_XY(url); % with output argument to just get the URL

XY  =mf_GMLL2XY(url.path.LATLON,url.center);

plot(XY(:,1),XY(:,2),'ro');

