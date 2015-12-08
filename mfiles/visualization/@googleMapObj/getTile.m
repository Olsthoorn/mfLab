function getTile(o,dLat,dLon,Lat,Lon)
%TODO TO 140110
% Let us retrieve an image from Google Maps of given resolution and size.
% The API:  Image = googleMapObj.getTile(dLat,dLon,Lat,Lon,zoomLevel)

% procedure, taking into account that google will at most issue an image of
% 640 by 640 pixels aound the point of interest at the requested zoom
% level.

% find out how many tiles we need around the point and which tiles to
% aquire
% Get title where center is in

% Get size of tile where center is in
% Compute the tiles where the four corners are in.
% Compute the coordinates inside the tiles.
latCorner = Lat + dLat * [-0.5; 0.5];
lonCorner = Lon + dLon + [-0.5, 0.5];

GMpoint = o(latCorner([1 1 2 2]),lonCorner([1 2 2 1]));

tiles = getTiles(GMpoint,zoomLevel);

% Joint the tiles
% if tiles.ix = tiles.ix, ...
%         if tiles.iy = tils.ix,...
%         end
% end

% TO140109 % too tired to continue

% Cout out the desired portion

