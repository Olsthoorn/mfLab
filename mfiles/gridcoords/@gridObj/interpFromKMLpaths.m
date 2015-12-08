function Layer = interpFromKMLpaths(o,Prefix)
%GRIDOBJ/interpFromKMLpaths(o,Prefix) -- spatially interpolate points
% based on isolines that have been digitized in Google Earth.
%
% A set of KML files containing Google Earth paths should be present
% in the directory and the ones to be used in the interpolation should
% all start with the given prefix followed by a nuber that corresponds
% with the contour line elevation. These file may for instance be called
% Depth00.kml, Depth05.kml Depth10.kml, Depth05b.kml etc. As long as they
% start with the given prefix and as long as that prefix is followed by a
% number that gives the thickness or depth of whatever. That number will be
% taken the elevation of the contourline contained in the file.
% The digitized points will be interpolated into a denser line with point
% distance set uniformly to 100 m.
%
% Then all the so-obtained points will be used to generate a spatial
% interpolator that will be applied to the grid to obtain an interpolated
% field whos coordinates correspond with the gridObj.
%
% Points outside the hull around the contour data will be NaN and
% correspond to areas that where no data area available.
%
% TO 130922

[P,Prefix] = fileparts(Prefix);

if isempty(Prefix), Prefix=P; P=''; end

if isempty(P)
    d  = dir([Prefix '*.kml']);
else
    d = dir([P filesep Prefix '*.kml']);
end

ds = 100;  % distance in degrees

% get the points
if isempty(P)
    Points = kmlFolder(d.name);
else
    Points = kmlFolder([P filesep d.name]);
end
for i=numel(Points):-1:1
    I =  (diff(Points(i).E)==0 & diff(Points(i).N)==0);
    if  any(I)
        Points(i).E(I) = [];
        Points(i).N(I) = [];
        Points(i).X(I) = [];
        Points(i).Y(I) = [];
    end
    try
        Points(i).Z = ones(size(Points(i).X))*sscanf(Points(i).name,[Prefix '%d']);
    catch  ME
        error([ME.message,'\nCan''t find thickness in contour <<%s>>'],Points(i).name);
        %Points(i)=[];
        %continue;
    end
    
    % compute length along line
    s1 = [0;  cumsum(sqrt(diff(Points(i).X).^2+diff(Points(i).Y).^2))];
    
    % generate equidistant points at distance ds
    s2 = [0:ds:s1(end) s1(end)]';
    
    Points(i).e = interp1(s1,Points(i).E,s2);
    Points(i).n = interp1(s1,Points(i).N,s2);
    [Points(i).x,Points(i).y] = wgs2utm(Points(i).n,Points(i).e);
    Points(i).z = Points(i).Z(1) *ones(size(Points(i).x));
end    

[x,y] = wgs2utm(vertcat(Points.n),vertcat(Points.e));

F = TriScatteredInterp(x,y,vertcat(Points.z));

% interpolate and make sure the MINDZ is satisfied

Layer = F(o.Xm,o.Ym);

Layer(isnan(Layer) | Layer<o.MINDZ) = o.MINDZ;


