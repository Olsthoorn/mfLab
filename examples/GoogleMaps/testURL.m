% TESTURL:  Requests a map directly from Google Maps into Matlab figure
%
% USAGE:
%   testURL
%
% DETAILS:
%   The functions can be found in mflab/mfiles/visualization
%
% TO 110506

close all;

%% Coordainats of an island from Dubai World
locations={...
[ 25.253159, 55.178000], 'Dubai-world island';
[ 52.372821,  4.893667], 'Dam Monument Amsterdam';
[ 68.962970, 33.089563], 'Murmansk, Russia';
[-33.856857,151.215192], 'Sydney';
[-54.795444,-68.232218], 'Ushuaia, Argentina';
[ 37.808810,-122.409803],'San Francisco, Fisherman''s Wharf';
[ 40.748524,-73.985676], 'New York Empire State building'
};

locs = cat(1,locations{:,1});

zoom    = 16;
pixels  = [512,512];
maptype = 'satellite';

clear GM;

for iLoc=size(locations,1):-1:1

    GM(iLoc) = googleMapObj(locs(iLoc,:),zoom,pixels,maptype,'png');
    
    GM(iLoc).name = locations{iLoc,2};

    GM(iLoc) = GM(iLoc).getImage();
    GM(iLoc).image([0 GM(iLoc).tile.Lx], [0 GM(iLoc).tile.Ly]);
    
    %% plot a reference grid
     xGr=0:250:1500; 
     yGr=0:250:1500;
     plotgrid(xGr,yGr,'w');
end

test_world
