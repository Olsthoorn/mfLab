% TEST_GM2PNG: EXAMPLE of how to use mf_GM2PNG to retrieve figure from Google Maps
%
% USAGE:
%    test_GM2PNG
%
% TO 110501 110514 121008

center='The+Palm+Jebel+Ali,Dubai,AE';    % as name
%center=[25.003395,54.986850];           % as coordinates

zoom    = 12;           % between 0 and 21
pixels  = [512 512]';   % each between 1 and 640
maptype = 'roadmap'; % one of

%% Add as many markers to the struct as you like
i=0;

i=i+1;
markers(i).color='blue';
markers(i).label='S';
markers(i).LATLON=[25.002696,54.988728];

i=i+1;
markers(i).color='green';
markers(i).label='G';
markers(i).LATLON=[ 25.030855, 54.971621];

i=i+1;
markers(i).color='red';
markers(i).label='C';
markers(i).LATLON=[ 24.982417,54.970744];

%% Add as many paths to your struct as you like
paths(1).color='blue';         % path color
paths(1).weight=5;             % path line thickness
paths(1).fillcolor='yellow';   % color to fill interior of path polygon

% actual path, as Lat Lon in decimal degrees
% you can use kmlpath('kmlpathfile.kml')  instead, where 'kmlpathfile.kml
% is the path clicked in Google Earth and saved as kml file.
paths(1).LATLON=[
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

GM = googleMapObj(center,zoom,pixels,maptype,'png',markers,paths);
GM.image;

