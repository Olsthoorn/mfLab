function [Lat,Lon]=utm2wgs(varargin)
%UTM2WGS converts the vectors of UTM coordinates into Lat/Lon vectors.
%
% Usage:
%   [Lat,Lon]=utm2wgs(xx,yy,utmzone)
%
% the UTM zone is the 6deg wide longitude zone counting eastward from te
% date line with a letter counting northward C to X. Notice that some zones
% are locally wider in the very north and south zones to incorporate a complete country.

% Inputs:
%    x       - UTM easting in meters
%    y       - UTM northing in meters
%    utmzone - UTM longitudinal zone
% Outputs:
%    Lat (WGS84 Latitude vector)  in decimal degrees:  ddd.dddddddd
%    Lon (WGS84 Longitude vector) in decimal degrees:  ddd.dddddddd
%
% Example:
%     x=[ 458731;  407653;  239027;  362850];
%     y=[4462881; 3171843; 4302285; 2772478];
%     utmzone=['30T'; '32T'; '01S'; '51R'];
%    [Lat,Lon]=utm2wgs(x,y,utmzone);
%       returns
% Lat =
%    40.3154
%    28.6705
%    38.8307
%    25.0618
% Lon =
%    -3.4857
%     8.0549
%  -180.0064
%   121.6403
%
% Source: DMA Technical Manual 8358.2, Fairfax, VA
%
% TO 110101

%% Argument checking

[xx,varargin] = getNext(varargin,'double',[]);
[yy,varargin] = getNext(varargin,'double',[]);
[utmzone,  ~] = getNext(varargin,{'char','cell'},[]);

if isempty(xx) ||isempty(yy) || isempty(utmzone)
    error(['%s: input should be xUTM,yUTM,utmzone(s), with utmzone(s) a string like ''30N''\n',...
        'an darray of such strings, or a cell array with such strings'],mfilename);
end

if ischar(utmzone)
    for i=size(utmzone,1):-1:1
        UTMZONE{i}=utmzone(i,:);
    end
else
    if iscell(utmzone)
        UTMZONE=utmzone;
    end
end

if numel(UTMZONE)>1
    j = numel(UTMZONE);
    for i=numel(xx):-1:j
        UTMZONE(i)=UTMZONE(j);
    end
end

if any(size(xx)~=size(yy))
   error('x, y must have same size');
end
c=size(UTMZONE{1},2);
if (c~=3)
   error('utmzone should be a vector of strings like "30T"');
end

%% Computing Lat/Lon coordinates for each input

for i=numel(xx):-1:1
    jj = min(numel(UTMZONE),i);
    if (UTMZONE{jj}(end)>'X' || UTMZONE{jj}(end)<'C')
        fprintf('utm2wgs: Warning you cannot use lowercase letters in UTM zone\n');
    end
    if UTMZONE{jj}(end)>'M'
        hemis='n';    % Northern hemisphere
    else
        hemis='s';    % Southern hemisphere
    end

    x=xx(i);
    y=yy(i);
    zone=str2double(UTMZONE{jj}(1:2));
    sa = 6378137.000000;                % semi-major axis of the Earth ellipsoid
    sb = 6356752.314245;                % semi-minor axis of the Earth ellipsoid
    e=(((sa^2)-(sb^2))^0.5)/sb;         % squared second eccentricity
    e2= e^2;
    c=sa^2/sb;
    X = x - 500000;
    if  hemis=='s'
        Y=y-10000000;
    else
        Y=y;
    end

    S=((zone*6)-183);
    lat=Y/(6366197.724*0.9996);
    v=(c/((1+(e2*(cos(lat))^2)))^0.5)*0.9996;
    a=X/v;
    a1=sin(2*lat);
    a2=a1*(cos(lat))^2;
    j2=lat+(a1/2);
    j4=((3*j2)+a2)/4;
    j6=((5*j4)+(a2*(cos(lat))^2))/3;
    alpha=(3/4)*e2;
    beta=(5/3)*alpha^2;
    gamma=(35/27)*alpha^3;
    Bm=0.9996*c*(lat-alpha*j2+beta*j4-gamma*j6);
    b=(Y-Bm)/v;
    Epsi=((e2*a^2)/2)*(cos(lat))^2;
    Eps=a*(1-(Epsi/3));
    nab=(b*(1-Epsi))+lat;
    senoheps=(exp(Eps)-exp(-Eps))/2;
    Delta=atan(senoheps/(cos(nab)));
    TaO=atan(cos(Delta)*tan(nab));
    longitude=(Delta*(180/pi))+S;
    latitude=(lat+(1+e2*(cos(lat)^2)-(3/2)*e2*sin(lat)*...
        cos(lat)*(TaO-lat))*(TaO-lat))*(180/pi);
    Lat(i)=latitude;
    Lon(i)=longitude;
end                 % For-loop end 
Lat=reshape(Lat,size(xx));
Lon=reshape(Lon,size(yy));