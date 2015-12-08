% example of how to retrieve a cross section from TNO REGIS


%getDinoXSec('mypath.kml');         % get XSection form dino loket along this Google Earth path
%getDinoXSec('WO grens B-NB.kml');  % get XSection form dino loket along this Google Earth path
getDinoXSec('AWD-WE.kml');          % get XSection form dino loket along this Google Earth path

%getRegisXS('mypath.kml'); % same but in the local Matlab web browser
%getRegisXS('WO grens B-NB.kml'); % same but in the local Matlab web browser

% I would like to retrieve the cross section direclty into a Matlab image
% however, without more details about how to request it, it seems
% difficult. Thanks to the documentation there is no such difficulty with
% respect to Google Maps figures, where the picture format is included in the request.
%  =====  [status,f]=urlwrite(URL, 'samples.html');
% In that case we can retrieve the figure

%A    = imread(URL);                 % get it directly into matlab
%info = imfinfo(URL);               % get associated image info
%image(xLim,yLim([2 1]),A); colormap(info.Colormap); % Plot it using the colormap from the info

% Here we can save the website. (right button, save as). We then obtain all
% files consituting the site, which includes the png files of the cross
% section and separately of its legend.
% TO 110514
